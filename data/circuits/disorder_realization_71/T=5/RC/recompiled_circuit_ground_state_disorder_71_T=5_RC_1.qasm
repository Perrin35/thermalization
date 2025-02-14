OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.355964) q[0];
sx q[0];
rz(-2.845919) q[0];
sx q[0];
rz(2.7276738) q[0];
rz(0.26894459) q[1];
sx q[1];
rz(-0.36965951) q[1];
sx q[1];
rz(-1.2764021) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3901509) q[0];
sx q[0];
rz(-0.12889114) q[0];
sx q[0];
rz(0.78123631) q[0];
rz(-pi) q[1];
rz(0.90246713) q[2];
sx q[2];
rz(-1.5215708) q[2];
sx q[2];
rz(0.38429457) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7404186) q[1];
sx q[1];
rz(-2.310862) q[1];
sx q[1];
rz(2.7526345) q[1];
rz(-1.7539361) q[3];
sx q[3];
rz(-2.1570342) q[3];
sx q[3];
rz(-0.30979118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74497491) q[2];
sx q[2];
rz(-1.2367542) q[2];
sx q[2];
rz(-1.2006987) q[2];
rz(-2.4602304) q[3];
sx q[3];
rz(-1.6187637) q[3];
sx q[3];
rz(2.7549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.50614) q[0];
sx q[0];
rz(-2.8400087) q[0];
sx q[0];
rz(0.26148456) q[0];
rz(-1.6188949) q[1];
sx q[1];
rz(-0.50917429) q[1];
sx q[1];
rz(-1.2759298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7000293) q[0];
sx q[0];
rz(-1.6671379) q[0];
sx q[0];
rz(-0.55192134) q[0];
x q[1];
rz(-0.65211703) q[2];
sx q[2];
rz(-2.0154833) q[2];
sx q[2];
rz(2.7114078) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2535227) q[1];
sx q[1];
rz(-1.2467398) q[1];
sx q[1];
rz(-0.56323207) q[1];
rz(-pi) q[2];
rz(-1.7959782) q[3];
sx q[3];
rz(-2.9888267) q[3];
sx q[3];
rz(-1.6880715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2168938) q[2];
sx q[2];
rz(-2.108722) q[2];
sx q[2];
rz(0.72803289) q[2];
rz(-1.7714455) q[3];
sx q[3];
rz(-2.5244505) q[3];
sx q[3];
rz(3.106015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034123357) q[0];
sx q[0];
rz(-2.0296445) q[0];
sx q[0];
rz(-2.4566101) q[0];
rz(2.0383539) q[1];
sx q[1];
rz(-2.5220242) q[1];
sx q[1];
rz(-1.0122976) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2489126) q[0];
sx q[0];
rz(-1.9143701) q[0];
sx q[0];
rz(-1.1569886) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9673789) q[2];
sx q[2];
rz(-1.763134) q[2];
sx q[2];
rz(-1.1735553) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.040739487) q[1];
sx q[1];
rz(-1.3636989) q[1];
sx q[1];
rz(-1.7000292) q[1];
x q[2];
rz(1.0676882) q[3];
sx q[3];
rz(-0.87260428) q[3];
sx q[3];
rz(-3.0029135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.181902) q[2];
sx q[2];
rz(-1.5519374) q[2];
sx q[2];
rz(0.8054136) q[2];
rz(2.3253333) q[3];
sx q[3];
rz(-2.5160242) q[3];
sx q[3];
rz(-0.094154112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020298088) q[0];
sx q[0];
rz(-2.2327406) q[0];
sx q[0];
rz(2.5138309) q[0];
rz(-0.10920564) q[1];
sx q[1];
rz(-0.86995482) q[1];
sx q[1];
rz(0.83762082) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9484884) q[0];
sx q[0];
rz(-2.645506) q[0];
sx q[0];
rz(-3.0015465) q[0];
x q[1];
rz(1.2455279) q[2];
sx q[2];
rz(-2.5211689) q[2];
sx q[2];
rz(1.1837981) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1091812) q[1];
sx q[1];
rz(-2.4285168) q[1];
sx q[1];
rz(0.20532691) q[1];
rz(0.67879642) q[3];
sx q[3];
rz(-0.46884109) q[3];
sx q[3];
rz(1.9285338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8236905) q[2];
sx q[2];
rz(-0.48419848) q[2];
sx q[2];
rz(-0.34453264) q[2];
rz(1.4065929) q[3];
sx q[3];
rz(-2.78648) q[3];
sx q[3];
rz(-3.0936892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59871167) q[0];
sx q[0];
rz(-0.6364091) q[0];
sx q[0];
rz(0.5058381) q[0];
rz(-2.089031) q[1];
sx q[1];
rz(-0.86741766) q[1];
sx q[1];
rz(0.94598407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.036412843) q[0];
sx q[0];
rz(-2.7650021) q[0];
sx q[0];
rz(2.3947477) q[0];
rz(-pi) q[1];
rz(0.35752488) q[2];
sx q[2];
rz(-1.4357743) q[2];
sx q[2];
rz(-2.6808878) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2767503) q[1];
sx q[1];
rz(-2.4227243) q[1];
sx q[1];
rz(2.5548773) q[1];
x q[2];
rz(0.068377062) q[3];
sx q[3];
rz(-1.7192327) q[3];
sx q[3];
rz(-2.1297586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.08789173) q[2];
sx q[2];
rz(-1.0975857) q[2];
sx q[2];
rz(1.5778479) q[2];
rz(2.0569233) q[3];
sx q[3];
rz(-2.2659149) q[3];
sx q[3];
rz(2.5904371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5797193) q[0];
sx q[0];
rz(-2.5358574) q[0];
sx q[0];
rz(2.9015923) q[0];
rz(2.2198246) q[1];
sx q[1];
rz(-1.0488989) q[1];
sx q[1];
rz(1.1704495) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83607996) q[0];
sx q[0];
rz(-1.5275947) q[0];
sx q[0];
rz(0.57010285) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7851849) q[2];
sx q[2];
rz(-1.6898167) q[2];
sx q[2];
rz(1.0072034) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45615087) q[1];
sx q[1];
rz(-0.69760453) q[1];
sx q[1];
rz(0.82932034) q[1];
x q[2];
rz(0.097436029) q[3];
sx q[3];
rz(-1.9698922) q[3];
sx q[3];
rz(-0.68828415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74448284) q[2];
sx q[2];
rz(-2.3220671) q[2];
sx q[2];
rz(-2.5605555) q[2];
rz(2.2033612) q[3];
sx q[3];
rz(-2.5751028) q[3];
sx q[3];
rz(0.70403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0051603) q[0];
sx q[0];
rz(-0.67661023) q[0];
sx q[0];
rz(-2.6410979) q[0];
rz(-0.23807921) q[1];
sx q[1];
rz(-1.6903189) q[1];
sx q[1];
rz(2.4949825) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71599267) q[0];
sx q[0];
rz(-1.5690205) q[0];
sx q[0];
rz(-3.1067757) q[0];
x q[1];
rz(3.1076381) q[2];
sx q[2];
rz(-1.6852549) q[2];
sx q[2];
rz(-1.7781228) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.38234455) q[1];
sx q[1];
rz(-2.0539762) q[1];
sx q[1];
rz(0.13570857) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.929333) q[3];
sx q[3];
rz(-1.6286116) q[3];
sx q[3];
rz(1.1392347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87767345) q[2];
sx q[2];
rz(-2.4223902) q[2];
sx q[2];
rz(0.70362299) q[2];
rz(-1.2879114) q[3];
sx q[3];
rz(-1.0602919) q[3];
sx q[3];
rz(2.6357486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8877761) q[0];
sx q[0];
rz(-2.0075338) q[0];
sx q[0];
rz(-1.8164841) q[0];
rz(2.3690986) q[1];
sx q[1];
rz(-1.1224116) q[1];
sx q[1];
rz(2.966029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52161509) q[0];
sx q[0];
rz(-0.78536638) q[0];
sx q[0];
rz(0.26544858) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2261691) q[2];
sx q[2];
rz(-2.2594514) q[2];
sx q[2];
rz(1.4530393) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9055134) q[1];
sx q[1];
rz(-1.4314974) q[1];
sx q[1];
rz(2.0361316) q[1];
rz(0.14388559) q[3];
sx q[3];
rz(-0.7215313) q[3];
sx q[3];
rz(0.084621457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.10302155) q[2];
sx q[2];
rz(-1.0744289) q[2];
sx q[2];
rz(-2.0477022) q[2];
rz(-1.966018) q[3];
sx q[3];
rz(-2.3837377) q[3];
sx q[3];
rz(0.28483835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71838251) q[0];
sx q[0];
rz(-3.0960313) q[0];
sx q[0];
rz(-1.3257931) q[0];
rz(-2.6153053) q[1];
sx q[1];
rz(-0.86295366) q[1];
sx q[1];
rz(0.78904271) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6279052) q[0];
sx q[0];
rz(-2.9531678) q[0];
sx q[0];
rz(2.1112676) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5297935) q[2];
sx q[2];
rz(-0.42610301) q[2];
sx q[2];
rz(0.95624051) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4261158) q[1];
sx q[1];
rz(-1.2271235) q[1];
sx q[1];
rz(2.2261376) q[1];
rz(0.6700683) q[3];
sx q[3];
rz(-0.89691554) q[3];
sx q[3];
rz(-1.64627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9290756) q[2];
sx q[2];
rz(-2.2079461) q[2];
sx q[2];
rz(-2.2361501) q[2];
rz(-0.91600156) q[3];
sx q[3];
rz(-1.6986877) q[3];
sx q[3];
rz(-2.0974832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4395897) q[0];
sx q[0];
rz(-2.3083394) q[0];
sx q[0];
rz(2.2138017) q[0];
rz(1.0935498) q[1];
sx q[1];
rz(-2.5826726) q[1];
sx q[1];
rz(-0.13339001) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0722712) q[0];
sx q[0];
rz(-1.7786553) q[0];
sx q[0];
rz(-2.8448807) q[0];
rz(-1.8850408) q[2];
sx q[2];
rz(-0.5945328) q[2];
sx q[2];
rz(-2.9192215) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9815326) q[1];
sx q[1];
rz(-1.0602177) q[1];
sx q[1];
rz(2.4909944) q[1];
rz(-pi) q[2];
x q[2];
rz(1.978558) q[3];
sx q[3];
rz(-2.6490903) q[3];
sx q[3];
rz(1.6655079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6371969) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(1.2633598) q[2];
rz(2.0718306) q[3];
sx q[3];
rz(-2.1061149) q[3];
sx q[3];
rz(0.076816471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99061154) q[0];
sx q[0];
rz(-2.4438416) q[0];
sx q[0];
rz(-1.7846815) q[0];
rz(-1.6928584) q[1];
sx q[1];
rz(-2.3285463) q[1];
sx q[1];
rz(2.9978233) q[1];
rz(-1.0634546) q[2];
sx q[2];
rz(-2.5567071) q[2];
sx q[2];
rz(-1.1699789) q[2];
rz(-1.9142022) q[3];
sx q[3];
rz(-1.8029984) q[3];
sx q[3];
rz(2.3912994) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
