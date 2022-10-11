function [loc] = genes_to_idxs(gene_list, genes, nmets)
    [~,loc] = ismember(gene_list,genes);
    loc = loc + nmets;
end
