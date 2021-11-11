select
'FT   CDS             complement('  || gene.pos_begin  || '..'  || gene.pos_end || ')',
'FT                   /locus_tag="' || gene.systematic_id ||'"',
'FT                   /gene="' || curated.name || '"',
'FT                   /product="' || curated.product || '"' ,
'FT                   /similarity="Similar to '|| similarto.organismo ||', '|| similarto.produto || ' (' || similarto.tamanho_subject || ' aa), e-value: '|| similarto.evalue ||', '|| similarto.percentual ||' id in '|| similarto.tamanho_query ||' aa"' ,
'FT                   /colour=' || decidecolor(gene.pseudogene, besthits.identity, besthits.size, gene.cds_size),
'FT                   /note="Local_subcelular(Surfg): ' || gene.local_subcelular || '"' ,
'FT                   /blastp_file="blastp/Cp19.embl.seq.0' || substring(gene.systematic_id, '....$') || '.out";',
getalldomain(gene.systematic_id),getsignal(gene.systematic_id),gettmh(gene.systematic_id)
from gene	LEFT OUTER JOIN (select distinct * from gene JOIN blasthits ON (gene.systematic_id = blasthits.query)
							 where blasthits.identity > 60 and blasthits.size > (gene.cds_size*0.60)
							) as besthits ON (gene.systematic_id = besthits.query)
			LEFT OUTER JOIN similarto  ON (gene.systematic_id = similarto.query)
			LEFT OUTER JOIN curated USING (subject)
where gene.orientation='-'
UNION
select
'FT   CDS             '  || gene.pos_begin  || '..'  || gene.pos_end,
'FT                   /locus_tag="' || gene.systematic_id ||'"',
'FT                   /gene="' || curated.name || '"',
'FT                   /product="' || curated.product || '"' ,
'FT                   /similarity="Similar to '|| similarto.organismo ||', '|| similarto.produto || ' (' || similarto.tamanho_subject || ' aa), e-value: '|| similarto.evalue ||', '|| similarto.percentual ||' id in '|| similarto.tamanho_query ||' aa"' ,
'FT                   /colour=' || decidecolor(gene.pseudogene, besthits.identity, besthits.size, gene.cds_size),
'FT                   /note="Local_subcelular(Surfg): ' || gene.local_subcelular || '"' ,
'FT                   /blastp_file="blastp/Cp19.embl.seq.0' || substring(gene.systematic_id, '....$') || '.out";',
getalldomain(gene.systematic_id),getsignal(gene.systematic_id),gettmh(gene.systematic_id)
from gene	LEFT OUTER JOIN (select distinct * from gene JOIN blasthits ON (gene.systematic_id = blasthits.query)
							 where blasthits.identity > 60 and blasthits.size > (gene.cds_size*0.60)
							) as besthits ON (gene.systematic_id = besthits.query)
			LEFT OUTER JOIN similarto  ON (gene.systematic_id = similarto.query)						
			LEFT OUTER JOIN curated USING (subject)
where gene.orientation='+'
