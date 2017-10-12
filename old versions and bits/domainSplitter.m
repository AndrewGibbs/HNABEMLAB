function splitPoints = domainSplitter( domain ) 
%splits domain up into subdomains, sort of
        splitPoints{1}=domain{1};
        splitPoints{2}=domain{2};
        if domain{1}(1)<domain{2}(2) && domain{1}(1)>domain{2}(1)
            splitPoints{2}=[splitPoints{2} domain{1}(1)];
        end
        if domain{1}(2)< domain{2}(2) && domain{1}(2)>domain{2}(1)
            splitPoints{2}=[splitPoints{2} domain{1}(2)];
        end

        if domain{2}(1)<domain{1}(2) && domain{2}(1)>domain{1}(1)
            splitPoints{1}=[splitPoints{1} domain{2}(1)];
        end
        if domain{2}(2)< domain{1}(2) && domain{2}(2)>domain{1}(1)
            splitPoints{1}=[splitPoints{1} domain{2}(2)];
        end

        splitPoints{1}=sort(union(domain{1}, splitPoints{1}));
        splitPoints{2}=sort(union(domain{2}, splitPoints{2}));

end

