OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.69230429) q[0];
sx q[0];
rz(-0.46625724) q[0];
sx q[0];
rz(1.847108) q[0];
rz(2.8757088) q[1];
sx q[1];
rz(2.9335913) q[1];
sx q[1];
rz(7.5367422) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82845518) q[0];
sx q[0];
rz(-1.0424409) q[0];
sx q[0];
rz(1.7067451) q[0];
rz(2.3135548) q[2];
sx q[2];
rz(-0.56664213) q[2];
sx q[2];
rz(0.41929454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47176018) q[1];
sx q[1];
rz(-1.3304632) q[1];
sx q[1];
rz(2.2369948) q[1];
x q[2];
rz(0.01913602) q[3];
sx q[3];
rz(-1.8950671) q[3];
sx q[3];
rz(-2.0024042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0365389) q[2];
sx q[2];
rz(-1.4194856) q[2];
sx q[2];
rz(1.6671906) q[2];
rz(0.27753943) q[3];
sx q[3];
rz(-1.8614635) q[3];
sx q[3];
rz(2.2138219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14811806) q[0];
sx q[0];
rz(-0.19522218) q[0];
sx q[0];
rz(-1.3268205) q[0];
rz(0.38816342) q[1];
sx q[1];
rz(-1.5048985) q[1];
sx q[1];
rz(-0.51663748) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9680133) q[0];
sx q[0];
rz(-1.359419) q[0];
sx q[0];
rz(-1.7239991) q[0];
x q[1];
rz(-2.0031163) q[2];
sx q[2];
rz(-0.45589635) q[2];
sx q[2];
rz(-1.0026463) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.109546) q[1];
sx q[1];
rz(-1.418772) q[1];
sx q[1];
rz(-0.8117453) q[1];
rz(2.2932986) q[3];
sx q[3];
rz(-2.775421) q[3];
sx q[3];
rz(-0.90400212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8119729) q[2];
sx q[2];
rz(-0.5054349) q[2];
sx q[2];
rz(0.86414117) q[2];
rz(-2.4498074) q[3];
sx q[3];
rz(-1.1312048) q[3];
sx q[3];
rz(-2.2085021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.316204) q[0];
sx q[0];
rz(-0.80556691) q[0];
sx q[0];
rz(1.4027931) q[0];
rz(-2.5750776) q[1];
sx q[1];
rz(-1.6367876) q[1];
sx q[1];
rz(-1.8623955) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4124547) q[0];
sx q[0];
rz(-1.8304145) q[0];
sx q[0];
rz(0.0065064001) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6635936) q[2];
sx q[2];
rz(-2.1076116) q[2];
sx q[2];
rz(-1.2144685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1058987) q[1];
sx q[1];
rz(-1.1269242) q[1];
sx q[1];
rz(-3.1015293) q[1];
rz(-1.1899105) q[3];
sx q[3];
rz(-2.1614868) q[3];
sx q[3];
rz(-2.4209972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0731571) q[2];
sx q[2];
rz(-2.7918039) q[2];
sx q[2];
rz(-1.5901828) q[2];
rz(2.6326211) q[3];
sx q[3];
rz(-0.83566982) q[3];
sx q[3];
rz(-1.6905748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.029025404) q[0];
sx q[0];
rz(-1.4968766) q[0];
sx q[0];
rz(-0.56370869) q[0];
rz(1.4504704) q[1];
sx q[1];
rz(-1.9861954) q[1];
sx q[1];
rz(0.82121003) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7696636) q[0];
sx q[0];
rz(-2.4255776) q[0];
sx q[0];
rz(0.016543702) q[0];
rz(-pi) q[1];
rz(-0.74521733) q[2];
sx q[2];
rz(-2.7181667) q[2];
sx q[2];
rz(-0.40027789) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5333818) q[1];
sx q[1];
rz(-2.1898263) q[1];
sx q[1];
rz(-2.4929895) q[1];
rz(-1.6333556) q[3];
sx q[3];
rz(-0.54397115) q[3];
sx q[3];
rz(2.6450752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.700909) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(1.7636501) q[2];
rz(-2.9772229) q[3];
sx q[3];
rz(-1.5956968) q[3];
sx q[3];
rz(0.36851287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57581562) q[0];
sx q[0];
rz(-1.6933279) q[0];
sx q[0];
rz(-0.50786316) q[0];
rz(2.2841618) q[1];
sx q[1];
rz(-1.4518167) q[1];
sx q[1];
rz(-0.47971496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0282766) q[0];
sx q[0];
rz(-1.2478267) q[0];
sx q[0];
rz(-0.66923176) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4998593) q[2];
sx q[2];
rz(-1.1878769) q[2];
sx q[2];
rz(-1.6429344) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8138201) q[1];
sx q[1];
rz(-1.346312) q[1];
sx q[1];
rz(2.1888365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8502838) q[3];
sx q[3];
rz(-1.5429792) q[3];
sx q[3];
rz(1.840855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5053284) q[2];
sx q[2];
rz(-2.0365066) q[2];
sx q[2];
rz(2.9218033) q[2];
rz(0.2291186) q[3];
sx q[3];
rz(-0.90355211) q[3];
sx q[3];
rz(2.1598099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55766469) q[0];
sx q[0];
rz(-0.54323498) q[0];
sx q[0];
rz(-1.1274717) q[0];
rz(0.30578956) q[1];
sx q[1];
rz(-1.9499648) q[1];
sx q[1];
rz(-2.9514899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8084304) q[0];
sx q[0];
rz(-1.8920533) q[0];
sx q[0];
rz(2.0684469) q[0];
x q[1];
rz(2.3476794) q[2];
sx q[2];
rz(-1.8474692) q[2];
sx q[2];
rz(-0.45396462) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87049228) q[1];
sx q[1];
rz(-1.2855347) q[1];
sx q[1];
rz(-2.3149957) q[1];
x q[2];
rz(0.28921952) q[3];
sx q[3];
rz(-2.8790124) q[3];
sx q[3];
rz(-2.5724831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1216792) q[2];
sx q[2];
rz(-1.4608258) q[2];
sx q[2];
rz(-1.100568) q[2];
rz(-1.3044283) q[3];
sx q[3];
rz(-1.5521939) q[3];
sx q[3];
rz(-1.7884375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8562427) q[0];
sx q[0];
rz(-0.42735639) q[0];
sx q[0];
rz(0.57619488) q[0];
rz(-2.088749) q[1];
sx q[1];
rz(-2.5710227) q[1];
sx q[1];
rz(2.0821234) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16082009) q[0];
sx q[0];
rz(-2.4376025) q[0];
sx q[0];
rz(-0.19613786) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9693107) q[2];
sx q[2];
rz(-1.5969807) q[2];
sx q[2];
rz(0.46253935) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6061343) q[1];
sx q[1];
rz(-2.6555057) q[1];
sx q[1];
rz(0.41907678) q[1];
rz(2.2031111) q[3];
sx q[3];
rz(-0.96358591) q[3];
sx q[3];
rz(-3.0724276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8656371) q[2];
sx q[2];
rz(-2.553678) q[2];
sx q[2];
rz(-1.0749344) q[2];
rz(2.2771207) q[3];
sx q[3];
rz(-1.266022) q[3];
sx q[3];
rz(1.5786494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39701715) q[0];
sx q[0];
rz(-0.96857849) q[0];
sx q[0];
rz(0.50719914) q[0];
rz(-1.1811258) q[1];
sx q[1];
rz(-1.6543417) q[1];
sx q[1];
rz(-2.2307253) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024689704) q[0];
sx q[0];
rz(-2.0729613) q[0];
sx q[0];
rz(2.2884918) q[0];
rz(-1.8751926) q[2];
sx q[2];
rz(-0.60303973) q[2];
sx q[2];
rz(1.837544) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6421318) q[1];
sx q[1];
rz(-1.5822463) q[1];
sx q[1];
rz(-1.7246555) q[1];
rz(-1.0491381) q[3];
sx q[3];
rz(-2.6781265) q[3];
sx q[3];
rz(-1.3090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38535038) q[2];
sx q[2];
rz(-1.4771947) q[2];
sx q[2];
rz(-0.64794668) q[2];
rz(-3.044965) q[3];
sx q[3];
rz(-2.1164618) q[3];
sx q[3];
rz(-2.1881762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18428093) q[0];
sx q[0];
rz(-2.1537557) q[0];
sx q[0];
rz(-1.8010358) q[0];
rz(2.6181009) q[1];
sx q[1];
rz(-1.610787) q[1];
sx q[1];
rz(-1.9370105) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0055375) q[0];
sx q[0];
rz(-2.3135472) q[0];
sx q[0];
rz(2.8715747) q[0];
rz(-pi) q[1];
rz(2.2080293) q[2];
sx q[2];
rz(-2.8945403) q[2];
sx q[2];
rz(-2.8480094) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5960418) q[1];
sx q[1];
rz(-1.4117472) q[1];
sx q[1];
rz(-0.56175128) q[1];
rz(-pi) q[2];
rz(-2.0818578) q[3];
sx q[3];
rz(-1.2132702) q[3];
sx q[3];
rz(-0.75393576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.039310731) q[2];
sx q[2];
rz(-1.0264341) q[2];
sx q[2];
rz(2.8459876) q[2];
rz(0.11232703) q[3];
sx q[3];
rz(-1.5887518) q[3];
sx q[3];
rz(2.2789392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91117793) q[0];
sx q[0];
rz(-0.4998315) q[0];
sx q[0];
rz(-2.3590132) q[0];
rz(-0.29620194) q[1];
sx q[1];
rz(-1.240851) q[1];
sx q[1];
rz(-2.9414419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91630477) q[0];
sx q[0];
rz(-0.80919832) q[0];
sx q[0];
rz(-1.7960219) q[0];
x q[1];
rz(-2.1364983) q[2];
sx q[2];
rz(-1.3315233) q[2];
sx q[2];
rz(-2.4636961) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0072277) q[1];
sx q[1];
rz(-2.0910462) q[1];
sx q[1];
rz(0.84905973) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8067935) q[3];
sx q[3];
rz(-2.349268) q[3];
sx q[3];
rz(2.5522422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.062181648) q[2];
sx q[2];
rz(-2.9997885) q[2];
sx q[2];
rz(-1.0395435) q[2];
rz(-0.027035106) q[3];
sx q[3];
rz(-0.98374933) q[3];
sx q[3];
rz(2.1496617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95299245) q[0];
sx q[0];
rz(-2.47692) q[0];
sx q[0];
rz(1.1823786) q[0];
rz(1.4324808) q[1];
sx q[1];
rz(-1.2624337) q[1];
sx q[1];
rz(1.591325) q[1];
rz(-1.6867137) q[2];
sx q[2];
rz(-1.7210261) q[2];
sx q[2];
rz(-2.6287519) q[2];
rz(0.19743528) q[3];
sx q[3];
rz(-2.6876269) q[3];
sx q[3];
rz(-0.85314565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
