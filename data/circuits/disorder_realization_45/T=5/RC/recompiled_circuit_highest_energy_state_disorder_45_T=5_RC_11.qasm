OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.6765321) q[0];
sx q[0];
rz(3.3829843) q[0];
sx q[0];
rz(9.2925185) q[0];
rz(-0.013068696) q[1];
sx q[1];
rz(2.4234534) q[1];
sx q[1];
rz(9.443774) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.976124) q[0];
sx q[0];
rz(-0.69688334) q[0];
sx q[0];
rz(-1.9096309) q[0];
rz(1.3842409) q[2];
sx q[2];
rz(-1.4779519) q[2];
sx q[2];
rz(1.873119) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.63848313) q[1];
sx q[1];
rz(-1.3967525) q[1];
sx q[1];
rz(0.32764224) q[1];
x q[2];
rz(-2.9844445) q[3];
sx q[3];
rz(-1.3408608) q[3];
sx q[3];
rz(0.63369753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12048177) q[2];
sx q[2];
rz(-2.3349473) q[2];
sx q[2];
rz(-0.79180229) q[2];
rz(2.0558489) q[3];
sx q[3];
rz(-1.2191685) q[3];
sx q[3];
rz(-2.048548) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2708112) q[0];
sx q[0];
rz(-2.9830611) q[0];
sx q[0];
rz(0.32919163) q[0];
rz(-1.5022494) q[1];
sx q[1];
rz(-0.90198016) q[1];
sx q[1];
rz(0.42218581) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3361583) q[0];
sx q[0];
rz(-1.0350739) q[0];
sx q[0];
rz(2.0690919) q[0];
rz(-pi) q[1];
rz(1.3866502) q[2];
sx q[2];
rz(-1.0304385) q[2];
sx q[2];
rz(-0.042405142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77743769) q[1];
sx q[1];
rz(-1.5244487) q[1];
sx q[1];
rz(2.8175161) q[1];
rz(-1.4677736) q[3];
sx q[3];
rz(-2.3405082) q[3];
sx q[3];
rz(-0.57263206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4210356) q[2];
sx q[2];
rz(-1.1853508) q[2];
sx q[2];
rz(2.0015008) q[2];
rz(0.0044048443) q[3];
sx q[3];
rz(-1.5811698) q[3];
sx q[3];
rz(0.13979039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7797101) q[0];
sx q[0];
rz(-0.10040586) q[0];
sx q[0];
rz(-2.7959339) q[0];
rz(-1.0935621) q[1];
sx q[1];
rz(-2.3193017) q[1];
sx q[1];
rz(2.0872769) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38906583) q[0];
sx q[0];
rz(-1.9574165) q[0];
sx q[0];
rz(1.9324612) q[0];
rz(-pi) q[1];
rz(-3.0385706) q[2];
sx q[2];
rz(-0.81562519) q[2];
sx q[2];
rz(1.6896923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7400764) q[1];
sx q[1];
rz(-1.1187727) q[1];
sx q[1];
rz(-1.4999309) q[1];
rz(-1.6003688) q[3];
sx q[3];
rz(-0.66638744) q[3];
sx q[3];
rz(2.8886384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.67636079) q[2];
sx q[2];
rz(-2.5939442) q[2];
sx q[2];
rz(-0.54222995) q[2];
rz(2.9690361) q[3];
sx q[3];
rz(-1.4901284) q[3];
sx q[3];
rz(-1.1057314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8050103) q[0];
sx q[0];
rz(-1.7886826) q[0];
sx q[0];
rz(1.2611058) q[0];
rz(-2.005596) q[1];
sx q[1];
rz(-2.0295862) q[1];
sx q[1];
rz(-3.0109308) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9016584) q[0];
sx q[0];
rz(-1.3011509) q[0];
sx q[0];
rz(-1.1872227) q[0];
rz(-pi) q[1];
rz(-2.6351852) q[2];
sx q[2];
rz(-1.7676438) q[2];
sx q[2];
rz(-1.0723237) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7893204) q[1];
sx q[1];
rz(-2.6323458) q[1];
sx q[1];
rz(1.5058084) q[1];
x q[2];
rz(-1.9404802) q[3];
sx q[3];
rz(-2.2013328) q[3];
sx q[3];
rz(-0.66555221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7589492) q[2];
sx q[2];
rz(-3.068277) q[2];
sx q[2];
rz(2.9200413) q[2];
rz(-0.70212901) q[3];
sx q[3];
rz(-1.5789072) q[3];
sx q[3];
rz(2.4041972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5704983) q[0];
sx q[0];
rz(-1.2449188) q[0];
sx q[0];
rz(0.13394295) q[0];
rz(-0.36930034) q[1];
sx q[1];
rz(-1.1993273) q[1];
sx q[1];
rz(1.7514924) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62890118) q[0];
sx q[0];
rz(-0.57399625) q[0];
sx q[0];
rz(2.0341134) q[0];
rz(-pi) q[1];
rz(-0.26525396) q[2];
sx q[2];
rz(-1.8657639) q[2];
sx q[2];
rz(-0.8910999) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.83727969) q[1];
sx q[1];
rz(-1.8472478) q[1];
sx q[1];
rz(2.5778092) q[1];
rz(-pi) q[2];
rz(-0.36365328) q[3];
sx q[3];
rz(-1.8580274) q[3];
sx q[3];
rz(1.938323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.2964581) q[2];
sx q[2];
rz(-2.239581) q[2];
sx q[2];
rz(1.5560537) q[2];
rz(-2.8499917) q[3];
sx q[3];
rz(-0.4762989) q[3];
sx q[3];
rz(2.8628023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45048243) q[0];
sx q[0];
rz(-2.1463558) q[0];
sx q[0];
rz(2.6158748) q[0];
rz(-2.763343) q[1];
sx q[1];
rz(-1.6098166) q[1];
sx q[1];
rz(-0.27413109) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4286684) q[0];
sx q[0];
rz(-0.98085058) q[0];
sx q[0];
rz(-1.601786) q[0];
rz(-pi) q[1];
rz(-2.3304105) q[2];
sx q[2];
rz(-0.59181442) q[2];
sx q[2];
rz(2.9331911) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25933274) q[1];
sx q[1];
rz(-1.3162398) q[1];
sx q[1];
rz(-1.987129) q[1];
rz(-pi) q[2];
rz(1.8242307) q[3];
sx q[3];
rz(-1.4259286) q[3];
sx q[3];
rz(1.0235909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5963001) q[2];
sx q[2];
rz(-1.2240852) q[2];
sx q[2];
rz(-2.8833585) q[2];
rz(1.5457414) q[3];
sx q[3];
rz(-1.7786547) q[3];
sx q[3];
rz(-2.7365007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79893583) q[0];
sx q[0];
rz(-2.6234143) q[0];
sx q[0];
rz(0.10547353) q[0];
rz(0.27490973) q[1];
sx q[1];
rz(-1.7173488) q[1];
sx q[1];
rz(-0.28604937) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3006004) q[0];
sx q[0];
rz(-1.6874041) q[0];
sx q[0];
rz(-2.9564315) q[0];
rz(-pi) q[1];
rz(2.3069165) q[2];
sx q[2];
rz(-1.3501985) q[2];
sx q[2];
rz(-2.5479864) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0487236) q[1];
sx q[1];
rz(-1.9561617) q[1];
sx q[1];
rz(3.0720622) q[1];
rz(-pi) q[2];
rz(2.4838832) q[3];
sx q[3];
rz(-2.4793787) q[3];
sx q[3];
rz(2.7194104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99453551) q[2];
sx q[2];
rz(-2.4592082) q[2];
sx q[2];
rz(0.4500173) q[2];
rz(0.58722812) q[3];
sx q[3];
rz(-1.3415895) q[3];
sx q[3];
rz(-1.2500866) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29379544) q[0];
sx q[0];
rz(-1.4601409) q[0];
sx q[0];
rz(-1.9548804) q[0];
rz(-0.31934357) q[1];
sx q[1];
rz(-1.0217383) q[1];
sx q[1];
rz(-0.26225463) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6410927) q[0];
sx q[0];
rz(-1.3670792) q[0];
sx q[0];
rz(0.55669703) q[0];
rz(-pi) q[1];
rz(-3.0056001) q[2];
sx q[2];
rz(-1.1609224) q[2];
sx q[2];
rz(0.52025822) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4353119) q[1];
sx q[1];
rz(-1.6526994) q[1];
sx q[1];
rz(1.0274853) q[1];
x q[2];
rz(-2.3877034) q[3];
sx q[3];
rz(-0.58710945) q[3];
sx q[3];
rz(-2.052149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.35683262) q[2];
sx q[2];
rz(-0.46669745) q[2];
sx q[2];
rz(-0.32877767) q[2];
rz(-1.6262936) q[3];
sx q[3];
rz(-1.3582151) q[3];
sx q[3];
rz(-0.40858194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(2.3560155) q[0];
sx q[0];
rz(-2.6308036) q[0];
sx q[0];
rz(-0.27895862) q[0];
rz(-0.51697671) q[1];
sx q[1];
rz(-2.7644283) q[1];
sx q[1];
rz(1.4252211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4734209) q[0];
sx q[0];
rz(-1.6145339) q[0];
sx q[0];
rz(-0.62181353) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4819229) q[2];
sx q[2];
rz(-2.2222509) q[2];
sx q[2];
rz(1.4209117) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.60630262) q[1];
sx q[1];
rz(-0.53724506) q[1];
sx q[1];
rz(-2.7657038) q[1];
x q[2];
rz(2.8513059) q[3];
sx q[3];
rz(-2.2072574) q[3];
sx q[3];
rz(-0.94062128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4031389) q[2];
sx q[2];
rz(-2.8454915) q[2];
sx q[2];
rz(-0.1599172) q[2];
rz(2.3219409) q[3];
sx q[3];
rz(-1.22217) q[3];
sx q[3];
rz(2.405449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4495471) q[0];
sx q[0];
rz(-1.8501546) q[0];
sx q[0];
rz(-0.29577574) q[0];
rz(-2.6462818) q[1];
sx q[1];
rz(-2.7073529) q[1];
sx q[1];
rz(-1.6326509) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2343649) q[0];
sx q[0];
rz(-1.1746658) q[0];
sx q[0];
rz(-2.7559126) q[0];
x q[1];
rz(3.0894439) q[2];
sx q[2];
rz(-1.9226769) q[2];
sx q[2];
rz(3.1089442) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1965643) q[1];
sx q[1];
rz(-1.6383645) q[1];
sx q[1];
rz(-3.0516207) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6466683) q[3];
sx q[3];
rz(-1.2240181) q[3];
sx q[3];
rz(2.1174255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8136924) q[2];
sx q[2];
rz(-1.6715965) q[2];
sx q[2];
rz(-0.63391614) q[2];
rz(3.0564803) q[3];
sx q[3];
rz(-1.2462933) q[3];
sx q[3];
rz(-2.3688721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2278628) q[0];
sx q[0];
rz(-1.5999404) q[0];
sx q[0];
rz(3.0318442) q[0];
rz(1.2211424) q[1];
sx q[1];
rz(-2.2962062) q[1];
sx q[1];
rz(-2.5795945) q[1];
rz(-0.22460266) q[2];
sx q[2];
rz(-1.9997659) q[2];
sx q[2];
rz(-1.275296) q[2];
rz(-1.3391277) q[3];
sx q[3];
rz(-1.5471519) q[3];
sx q[3];
rz(0.027761264) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
