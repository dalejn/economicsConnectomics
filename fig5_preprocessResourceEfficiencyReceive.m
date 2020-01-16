nreg=360;
nedge=64620;

cd('/data/jux/BBL/projects/ASLnetwork/data/faGlasser')
fa_network_files = dir('*FA_GlasserPNC.mat');
nfiles = length(fa_network_files);
nsub=nfiles;
fa_sq = zeros(nsub, nedge);

for k = 1:nfiles
    fa_net = load(fa_network_files(k).name);
    fa_net = fa_net.connectivity;
    fa_net = fa_net - diag(diag(fa_net));
    fa_sq(k,:) = squareform(fa_net);
end

cd /data/jux/BBL/projects/ASLnetwork/scripts
for probability=[0.999,0.98,0.96,0.94,0.92,0.9,0.8,0.7,0.6,0.5,0.4, 0.3, 0.2, 0.1]
    disp(probability)
    resEff_array=zeros(nsub,nreg+1);

    for s=1:nsub
        A_fa=squareform(fa_sq(s,:));
        A = A_fa - diag(diag(A_fa));

        [comps,~] = get_components(A);
        A = A(comps==1,comps==1);
        islands = find(comps>1);

        Eres = resource_efficiency_wei(A,probability);
        Eres = nanmean(Eres,1);

        if length(islands)>0
            for k=1:length(islands)
                i = islands(k);
                if k>1
                    B = [B(1:i-1),nan,B(i:end)];
                else
                    B = [Eres(1:i-1),nan,Eres(i:end)];
                end
            end
            resEff_array(s,2:nreg+1) = B;
        else 
            resEff_array(s,2:nreg+1) = Eres;
        end
    %    N = size(A,1);
    %    GEres = mean(Eres(~eye(N)>0))

        subj_name = fa_network_files(s).name
        resEff_array(s,1) = str2num(strtok(subj_name, '_'));
    end
    currentProb = num2str(probability, 16);
    currentProb(currentProb=='.') = []
    outName = strcat('/data/jux/BBL/projects/ASLnetwork/scripts/zaixuRepro/data/resource_efficiency',currentProb,'_receive.txt')
    dlmwrite(outName,resEff_array,'delimiter',' ', 'precision', 10)
end
