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
rz(-0.7402339) q[0];
sx q[0];
rz(-1.482168) q[0];
sx q[0];
rz(2.8066714) q[0];
rz(-2.6236293) q[1];
sx q[1];
rz(-2.1393175) q[1];
sx q[1];
rz(-0.60751539) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7142732) q[0];
sx q[0];
rz(-1.8022984) q[0];
sx q[0];
rz(-2.0718859) q[0];
rz(-pi) q[1];
rz(-3.0246441) q[2];
sx q[2];
rz(-2.0772572) q[2];
sx q[2];
rz(-3.1053271) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0820056) q[1];
sx q[1];
rz(-1.730189) q[1];
sx q[1];
rz(-2.6578147) q[1];
rz(-2.0161765) q[3];
sx q[3];
rz(-1.9102671) q[3];
sx q[3];
rz(-2.826863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.5241549) q[2];
sx q[2];
rz(-1.2158771) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(3.0607306) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(-1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2422159) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(2.2826165) q[0];
rz(-1.8513177) q[1];
sx q[1];
rz(-1.6764418) q[1];
sx q[1];
rz(1.4069517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1024723) q[0];
sx q[0];
rz(-2.6022096) q[0];
sx q[0];
rz(1.1155737) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036769899) q[2];
sx q[2];
rz(-1.2388133) q[2];
sx q[2];
rz(-2.9197502) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42900742) q[1];
sx q[1];
rz(-2.6204094) q[1];
sx q[1];
rz(-1.7527449) q[1];
rz(-pi) q[2];
rz(-2.5464393) q[3];
sx q[3];
rz(-0.66616026) q[3];
sx q[3];
rz(-0.43597886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0824288) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(-1.814369) q[2];
rz(2.0969157) q[3];
sx q[3];
rz(-0.71377126) q[3];
sx q[3];
rz(2.7248342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5663719) q[0];
sx q[0];
rz(-2.2307668) q[0];
sx q[0];
rz(-0.4253934) q[0];
rz(1.3771903) q[1];
sx q[1];
rz(-1.4981937) q[1];
sx q[1];
rz(1.4345217) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1578428) q[0];
sx q[0];
rz(-0.93261496) q[0];
sx q[0];
rz(1.4236091) q[0];
rz(0.038952053) q[2];
sx q[2];
rz(-2.4329081) q[2];
sx q[2];
rz(1.5295636) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9454398) q[1];
sx q[1];
rz(-1.1588133) q[1];
sx q[1];
rz(-0.71228551) q[1];
x q[2];
rz(2.8355153) q[3];
sx q[3];
rz(-0.97460213) q[3];
sx q[3];
rz(-1.5709189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7003358) q[2];
sx q[2];
rz(-1.2563027) q[2];
sx q[2];
rz(-1.2909935) q[2];
rz(-1.357632) q[3];
sx q[3];
rz(-1.6358401) q[3];
sx q[3];
rz(1.6443845) q[3];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62035471) q[0];
sx q[0];
rz(-1.0076032) q[0];
sx q[0];
rz(-2.3714016) q[0];
rz(-0.99984804) q[1];
sx q[1];
rz(-2.5412173) q[1];
sx q[1];
rz(1.8005449) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1091546) q[0];
sx q[0];
rz(-2.1103835) q[0];
sx q[0];
rz(1.3225609) q[0];
rz(-pi) q[1];
rz(-1.0275339) q[2];
sx q[2];
rz(-0.76342602) q[2];
sx q[2];
rz(-1.9486519) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0049147) q[1];
sx q[1];
rz(-1.305849) q[1];
sx q[1];
rz(0.30989225) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8982696) q[3];
sx q[3];
rz(-2.8411715) q[3];
sx q[3];
rz(-1.6302366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(2.3731903) q[2];
rz(0.10041222) q[3];
sx q[3];
rz(-1.5407591) q[3];
sx q[3];
rz(-1.8628619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95555821) q[0];
sx q[0];
rz(-2.8362507) q[0];
sx q[0];
rz(-0.22698639) q[0];
rz(1.3817894) q[1];
sx q[1];
rz(-2.558936) q[1];
sx q[1];
rz(-1.3963799) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5433301) q[0];
sx q[0];
rz(-2.6285183) q[0];
sx q[0];
rz(0.76900478) q[0];
rz(-2.7311205) q[2];
sx q[2];
rz(-1.7232401) q[2];
sx q[2];
rz(1.6688011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.760031) q[1];
sx q[1];
rz(-2.2027822) q[1];
sx q[1];
rz(-1.4820335) q[1];
rz(-pi) q[2];
rz(0.030524039) q[3];
sx q[3];
rz(-0.82477335) q[3];
sx q[3];
rz(1.8061226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.22623006) q[2];
sx q[2];
rz(-1.1504983) q[2];
sx q[2];
rz(3.0653595) q[2];
rz(-1.3876312) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(-1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6607894) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(-1.8060818) q[0];
rz(-2.238359) q[1];
sx q[1];
rz(-1.8025554) q[1];
sx q[1];
rz(-2.9551771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7714027) q[0];
sx q[0];
rz(-1.0433084) q[0];
sx q[0];
rz(-0.50748388) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30419402) q[2];
sx q[2];
rz(-0.59202164) q[2];
sx q[2];
rz(-1.5409868) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32115667) q[1];
sx q[1];
rz(-1.8119436) q[1];
sx q[1];
rz(0.92280573) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0669492) q[3];
sx q[3];
rz(-2.0130139) q[3];
sx q[3];
rz(2.7520455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.67941252) q[2];
sx q[2];
rz(-1.1932411) q[2];
sx q[2];
rz(-1.3235486) q[2];
rz(-1.9994252) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(-0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46397504) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(0.63419813) q[0];
rz(2.0967261) q[1];
sx q[1];
rz(-0.8909145) q[1];
sx q[1];
rz(-0.96010906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2506367) q[0];
sx q[0];
rz(-1.5433558) q[0];
sx q[0];
rz(-1.1720042) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5788001) q[2];
sx q[2];
rz(-1.8718613) q[2];
sx q[2];
rz(2.0640399) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83852488) q[1];
sx q[1];
rz(-2.6160598) q[1];
sx q[1];
rz(-1.5477033) q[1];
rz(2.3784901) q[3];
sx q[3];
rz(-2.1402485) q[3];
sx q[3];
rz(0.72009898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94379696) q[2];
sx q[2];
rz(-2.4891977) q[2];
sx q[2];
rz(1.5956399) q[2];
rz(1.4173077) q[3];
sx q[3];
rz(-2.3167819) q[3];
sx q[3];
rz(-1.6212757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5328131) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(-0.97484318) q[0];
rz(-3.0335562) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.1688165) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9038531) q[0];
sx q[0];
rz(-2.0444336) q[0];
sx q[0];
rz(-0.4717548) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5708099) q[2];
sx q[2];
rz(-2.332649) q[2];
sx q[2];
rz(-0.9521614) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7119693) q[1];
sx q[1];
rz(-1.6584466) q[1];
sx q[1];
rz(1.3369249) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6498472) q[3];
sx q[3];
rz(-1.5881268) q[3];
sx q[3];
rz(-0.31760707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.10890266) q[2];
sx q[2];
rz(-2.2638075) q[2];
sx q[2];
rz(2.4857793) q[2];
rz(-0.10410318) q[3];
sx q[3];
rz(-1.7963573) q[3];
sx q[3];
rz(-2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8769237) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(-1.4917829) q[0];
rz(-1.0393633) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(-0.46844354) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7130947) q[0];
sx q[0];
rz(-1.575043) q[0];
sx q[0];
rz(-2.9539786) q[0];
rz(2.6985964) q[2];
sx q[2];
rz(-1.3123242) q[2];
sx q[2];
rz(1.6719831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5891582) q[1];
sx q[1];
rz(-0.8048519) q[1];
sx q[1];
rz(-2.9167152) q[1];
rz(-pi) q[2];
rz(1.2896721) q[3];
sx q[3];
rz(-1.2996965) q[3];
sx q[3];
rz(-2.2106314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.046772) q[2];
sx q[2];
rz(-1.5609317) q[2];
sx q[2];
rz(-2.2136733) q[2];
rz(1.8038484) q[3];
sx q[3];
rz(-2.6313621) q[3];
sx q[3];
rz(-1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0874262) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(2.0637276) q[0];
rz(-0.36733356) q[1];
sx q[1];
rz(-1.3307738) q[1];
sx q[1];
rz(1.8574538) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3340942) q[0];
sx q[0];
rz(-1.1002514) q[0];
sx q[0];
rz(-0.28326359) q[0];
rz(1.7214891) q[2];
sx q[2];
rz(-0.65905276) q[2];
sx q[2];
rz(-2.6113752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2312647) q[1];
sx q[1];
rz(-2.961425) q[1];
sx q[1];
rz(2.641201) q[1];
rz(-0.73478847) q[3];
sx q[3];
rz(-1.7289203) q[3];
sx q[3];
rz(-1.8693592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.0014570634) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(0.4293116) q[2];
rz(0.15520994) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(-1.8163053) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24406381) q[0];
sx q[0];
rz(-1.5499935) q[0];
sx q[0];
rz(1.5503379) q[0];
rz(-0.86391972) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(-2.5979832) q[2];
sx q[2];
rz(-1.7805486) q[2];
sx q[2];
rz(-1.8303895) q[2];
rz(-0.33959099) q[3];
sx q[3];
rz(-0.97151269) q[3];
sx q[3];
rz(1.1480939) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
