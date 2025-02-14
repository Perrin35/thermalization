OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7856287) q[0];
sx q[0];
rz(-0.2956737) q[0];
sx q[0];
rz(0.41391882) q[0];
rz(3.4105372) q[1];
sx q[1];
rz(6.6528448) q[1];
sx q[1];
rz(11.289968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3901509) q[0];
sx q[0];
rz(-0.12889114) q[0];
sx q[0];
rz(-0.78123631) q[0];
rz(-pi) q[1];
rz(1.6501313) q[2];
sx q[2];
rz(-0.66986194) q[2];
sx q[2];
rz(-1.8928493) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40117404) q[1];
sx q[1];
rz(-2.310862) q[1];
sx q[1];
rz(0.38895815) q[1];
rz(-pi) q[2];
rz(-1.3876565) q[3];
sx q[3];
rz(-2.1570342) q[3];
sx q[3];
rz(-2.8318015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74497491) q[2];
sx q[2];
rz(-1.9048385) q[2];
sx q[2];
rz(-1.940894) q[2];
rz(-0.68136224) q[3];
sx q[3];
rz(-1.5228289) q[3];
sx q[3];
rz(2.7549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.50614) q[0];
sx q[0];
rz(-2.8400087) q[0];
sx q[0];
rz(-0.26148456) q[0];
rz(-1.5226978) q[1];
sx q[1];
rz(-2.6324184) q[1];
sx q[1];
rz(1.8656628) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0715213) q[0];
sx q[0];
rz(-1.0217279) q[0];
sx q[0];
rz(1.683805) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65211703) q[2];
sx q[2];
rz(-1.1261093) q[2];
sx q[2];
rz(2.7114078) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.88807) q[1];
sx q[1];
rz(-1.8948529) q[1];
sx q[1];
rz(0.56323207) q[1];
x q[2];
rz(0.034364465) q[3];
sx q[3];
rz(-1.4219163) q[3];
sx q[3];
rz(1.9158165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9246989) q[2];
sx q[2];
rz(-1.0328707) q[2];
sx q[2];
rz(2.4135598) q[2];
rz(1.7714455) q[3];
sx q[3];
rz(-2.5244505) q[3];
sx q[3];
rz(-3.106015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.034123357) q[0];
sx q[0];
rz(-2.0296445) q[0];
sx q[0];
rz(0.6849826) q[0];
rz(-1.1032387) q[1];
sx q[1];
rz(-0.61956844) q[1];
sx q[1];
rz(1.0122976) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023803614) q[0];
sx q[0];
rz(-2.6101681) q[0];
sx q[0];
rz(0.84367911) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2980897) q[2];
sx q[2];
rz(-2.8828104) q[2];
sx q[2];
rz(-2.7121787) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.040739487) q[1];
sx q[1];
rz(-1.7778938) q[1];
sx q[1];
rz(1.7000292) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0739044) q[3];
sx q[3];
rz(-2.2689884) q[3];
sx q[3];
rz(3.0029135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.95969069) q[2];
sx q[2];
rz(-1.5519374) q[2];
sx q[2];
rz(-2.336179) q[2];
rz(-2.3253333) q[3];
sx q[3];
rz(-2.5160242) q[3];
sx q[3];
rz(0.094154112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1212946) q[0];
sx q[0];
rz(-2.2327406) q[0];
sx q[0];
rz(-2.5138309) q[0];
rz(3.032387) q[1];
sx q[1];
rz(-0.86995482) q[1];
sx q[1];
rz(0.83762082) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9484884) q[0];
sx q[0];
rz(-2.645506) q[0];
sx q[0];
rz(3.0015465) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1659746) q[2];
sx q[2];
rz(-1.7576697) q[2];
sx q[2];
rz(0.65480168) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0324114) q[1];
sx q[1];
rz(-2.4285168) q[1];
sx q[1];
rz(-0.20532691) q[1];
rz(-pi) q[2];
rz(-1.8786976) q[3];
sx q[3];
rz(-1.9301722) q[3];
sx q[3];
rz(-1.1933768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8236905) q[2];
sx q[2];
rz(-0.48419848) q[2];
sx q[2];
rz(2.79706) q[2];
rz(-1.4065929) q[3];
sx q[3];
rz(-2.78648) q[3];
sx q[3];
rz(3.0936892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59871167) q[0];
sx q[0];
rz(-2.5051835) q[0];
sx q[0];
rz(-2.6357546) q[0];
rz(-1.0525616) q[1];
sx q[1];
rz(-0.86741766) q[1];
sx q[1];
rz(2.1956086) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3220468) q[0];
sx q[0];
rz(-1.2975386) q[0];
sx q[0];
rz(1.8332493) q[0];
rz(0.37028124) q[2];
sx q[2];
rz(-0.38114377) q[2];
sx q[2];
rz(-0.76424341) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.16984384) q[1];
sx q[1];
rz(-1.1976114) q[1];
sx q[1];
rz(-0.62974522) q[1];
rz(-pi) q[2];
rz(-1.4220174) q[3];
sx q[3];
rz(-1.5031723) q[3];
sx q[3];
rz(-0.54883445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0537009) q[2];
sx q[2];
rz(-1.0975857) q[2];
sx q[2];
rz(-1.5778479) q[2];
rz(-1.0846694) q[3];
sx q[3];
rz(-2.2659149) q[3];
sx q[3];
rz(2.5904371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56187335) q[0];
sx q[0];
rz(-2.5358574) q[0];
sx q[0];
rz(-0.24000034) q[0];
rz(2.2198246) q[1];
sx q[1];
rz(-1.0488989) q[1];
sx q[1];
rz(-1.9711432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3791948) q[0];
sx q[0];
rz(-1.0012915) q[0];
sx q[0];
rz(1.5194917) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7851849) q[2];
sx q[2];
rz(-1.6898167) q[2];
sx q[2];
rz(1.0072034) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7237548) q[1];
sx q[1];
rz(-1.0772632) q[1];
sx q[1];
rz(-2.626471) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3440873) q[3];
sx q[3];
rz(-2.7313958) q[3];
sx q[3];
rz(2.6997379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.74448284) q[2];
sx q[2];
rz(-0.8195256) q[2];
sx q[2];
rz(2.5605555) q[2];
rz(2.2033612) q[3];
sx q[3];
rz(-0.56648985) q[3];
sx q[3];
rz(-0.70403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0051603) q[0];
sx q[0];
rz(-0.67661023) q[0];
sx q[0];
rz(0.50049472) q[0];
rz(-0.23807921) q[1];
sx q[1];
rz(-1.6903189) q[1];
sx q[1];
rz(2.4949825) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4256) q[0];
sx q[0];
rz(-1.5690205) q[0];
sx q[0];
rz(0.034816936) q[0];
x q[1];
rz(-1.2836566) q[2];
sx q[2];
rz(-3.022225) q[2];
sx q[2];
rz(1.0743846) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.096488086) q[1];
sx q[1];
rz(-2.641171) q[1];
sx q[1];
rz(1.3184271) q[1];
x q[2];
rz(-0.061731652) q[3];
sx q[3];
rz(-1.2128856) q[3];
sx q[3];
rz(-0.45321143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.87767345) q[2];
sx q[2];
rz(-2.4223902) q[2];
sx q[2];
rz(-2.4379697) q[2];
rz(-1.2879114) q[3];
sx q[3];
rz(-2.0813007) q[3];
sx q[3];
rz(0.50584403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8877761) q[0];
sx q[0];
rz(-1.1340589) q[0];
sx q[0];
rz(1.8164841) q[0];
rz(2.3690986) q[1];
sx q[1];
rz(-1.1224116) q[1];
sx q[1];
rz(2.966029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52161509) q[0];
sx q[0];
rz(-2.3562263) q[0];
sx q[0];
rz(-2.8761441) q[0];
rz(-pi) q[1];
rz(-0.80412038) q[2];
sx q[2];
rz(-2.0607227) q[2];
sx q[2];
rz(-2.8049289) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0649291) q[1];
sx q[1];
rz(-2.6573219) q[1];
sx q[1];
rz(1.8736429) q[1];
rz(1.4453078) q[3];
sx q[3];
rz(-2.2832722) q[3];
sx q[3];
rz(-3.0355796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0385711) q[2];
sx q[2];
rz(-1.0744289) q[2];
sx q[2];
rz(-1.0938905) q[2];
rz(-1.1755747) q[3];
sx q[3];
rz(-2.3837377) q[3];
sx q[3];
rz(-0.28483835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4232101) q[0];
sx q[0];
rz(-3.0960313) q[0];
sx q[0];
rz(-1.3257931) q[0];
rz(-0.52628738) q[1];
sx q[1];
rz(-2.278639) q[1];
sx q[1];
rz(0.78904271) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1068971) q[0];
sx q[0];
rz(-1.409484) q[0];
sx q[0];
rz(0.09780305) q[0];
x q[1];
rz(-0.61179917) q[2];
sx q[2];
rz(-2.7154896) q[2];
sx q[2];
rz(2.1853521) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0328417) q[1];
sx q[1];
rz(-2.1819253) q[1];
sx q[1];
rz(-0.42400363) q[1];
rz(-pi) q[2];
rz(-2.2317722) q[3];
sx q[3];
rz(-2.2299521) q[3];
sx q[3];
rz(2.5499663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9290756) q[2];
sx q[2];
rz(-2.2079461) q[2];
sx q[2];
rz(-2.2361501) q[2];
rz(-2.2255911) q[3];
sx q[3];
rz(-1.4429049) q[3];
sx q[3];
rz(1.0441095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7020029) q[0];
sx q[0];
rz(-0.8332533) q[0];
sx q[0];
rz(-2.2138017) q[0];
rz(1.0935498) q[1];
sx q[1];
rz(-2.5826726) q[1];
sx q[1];
rz(-0.13339001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0722712) q[0];
sx q[0];
rz(-1.3629374) q[0];
sx q[0];
rz(0.29671191) q[0];
rz(1.2565518) q[2];
sx q[2];
rz(-2.5470599) q[2];
sx q[2];
rz(-0.22237117) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9815326) q[1];
sx q[1];
rz(-2.0813749) q[1];
sx q[1];
rz(-0.65059827) q[1];
rz(-0.20966704) q[3];
sx q[3];
rz(-2.0197967) q[3];
sx q[3];
rz(2.1213139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6371969) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(1.2633598) q[2];
rz(-1.0697621) q[3];
sx q[3];
rz(-2.1061149) q[3];
sx q[3];
rz(-3.0647762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99061154) q[0];
sx q[0];
rz(-0.69775109) q[0];
sx q[0];
rz(1.3569111) q[0];
rz(1.6928584) q[1];
sx q[1];
rz(-0.81304638) q[1];
sx q[1];
rz(-0.1437694) q[1];
rz(-1.0634546) q[2];
sx q[2];
rz(-2.5567071) q[2];
sx q[2];
rz(-1.1699789) q[2];
rz(-1.2273905) q[3];
sx q[3];
rz(-1.3385942) q[3];
sx q[3];
rz(-0.75029324) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
