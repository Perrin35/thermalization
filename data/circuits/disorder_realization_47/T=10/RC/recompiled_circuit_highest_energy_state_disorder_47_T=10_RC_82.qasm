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
rz(-2.0464309) q[0];
sx q[0];
rz(-0.14552966) q[0];
sx q[0];
rz(-2.6382883) q[0];
rz(-1.9637928) q[1];
sx q[1];
rz(4.7363321) q[1];
sx q[1];
rz(7.9793032) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2615525) q[0];
sx q[0];
rz(-1.6140811) q[0];
sx q[0];
rz(-0.39899428) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2022871) q[2];
sx q[2];
rz(-2.21647) q[2];
sx q[2];
rz(-1.1935357) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1481944) q[1];
sx q[1];
rz(-2.1752417) q[1];
sx q[1];
rz(-1.9733713) q[1];
rz(-pi) q[2];
rz(-1.5551994) q[3];
sx q[3];
rz(-1.5216148) q[3];
sx q[3];
rz(1.552207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4989) q[2];
sx q[2];
rz(-1.91012) q[2];
sx q[2];
rz(1.5042245) q[2];
rz(0.73689342) q[3];
sx q[3];
rz(-1.6156018) q[3];
sx q[3];
rz(-0.98114291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7351643) q[0];
sx q[0];
rz(-0.46877113) q[0];
sx q[0];
rz(0.51097393) q[0];
rz(1.3276395) q[1];
sx q[1];
rz(-1.4013441) q[1];
sx q[1];
rz(-0.32646349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8998979) q[0];
sx q[0];
rz(-2.7621671) q[0];
sx q[0];
rz(2.3510758) q[0];
x q[1];
rz(1.3148221) q[2];
sx q[2];
rz(-0.69967383) q[2];
sx q[2];
rz(0.2521317) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6447811) q[1];
sx q[1];
rz(-2.4678951) q[1];
sx q[1];
rz(2.069887) q[1];
rz(-0.28387733) q[3];
sx q[3];
rz(-1.9440158) q[3];
sx q[3];
rz(-2.9220327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1702801) q[2];
sx q[2];
rz(-2.9176517) q[2];
sx q[2];
rz(-1.2962606) q[2];
rz(2.7235459) q[3];
sx q[3];
rz(-0.97801912) q[3];
sx q[3];
rz(-0.70438284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.068785) q[0];
sx q[0];
rz(-0.3827706) q[0];
sx q[0];
rz(-2.5991154) q[0];
rz(-1.3756649) q[1];
sx q[1];
rz(-2.723697) q[1];
sx q[1];
rz(3.1105522) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.870793) q[0];
sx q[0];
rz(-2.7794771) q[0];
sx q[0];
rz(2.5848081) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.228187) q[2];
sx q[2];
rz(-1.3904461) q[2];
sx q[2];
rz(-3.0424398) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6185304) q[1];
sx q[1];
rz(-1.8162604) q[1];
sx q[1];
rz(-2.1851468) q[1];
x q[2];
rz(-2.5865745) q[3];
sx q[3];
rz(-1.5848098) q[3];
sx q[3];
rz(2.009258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5151908) q[2];
sx q[2];
rz(-2.4050737) q[2];
sx q[2];
rz(0.82702965) q[2];
rz(2.6229897) q[3];
sx q[3];
rz(-1.1122455) q[3];
sx q[3];
rz(1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9943635) q[0];
sx q[0];
rz(-1.8356859) q[0];
sx q[0];
rz(1.2116785) q[0];
rz(1.0481102) q[1];
sx q[1];
rz(-2.4217889) q[1];
sx q[1];
rz(-1.3406219) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7692684) q[0];
sx q[0];
rz(-2.4256673) q[0];
sx q[0];
rz(2.6628519) q[0];
rz(0.55669341) q[2];
sx q[2];
rz(-0.27678267) q[2];
sx q[2];
rz(-2.3900696) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2681899) q[1];
sx q[1];
rz(-2.4854069) q[1];
sx q[1];
rz(-0.49439685) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31130917) q[3];
sx q[3];
rz(-1.177586) q[3];
sx q[3];
rz(1.9462727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2010605) q[2];
sx q[2];
rz(-2.9643855) q[2];
sx q[2];
rz(1.1411257) q[2];
rz(-1.5308258) q[3];
sx q[3];
rz(-2.174236) q[3];
sx q[3];
rz(-0.78260666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.57946) q[0];
sx q[0];
rz(-1.9714404) q[0];
sx q[0];
rz(2.539047) q[0];
rz(0.58019477) q[1];
sx q[1];
rz(-1.6289026) q[1];
sx q[1];
rz(-0.7811195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21602042) q[0];
sx q[0];
rz(-1.4013774) q[0];
sx q[0];
rz(1.8623307) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65004543) q[2];
sx q[2];
rz(-2.263265) q[2];
sx q[2];
rz(-1.6672857) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8153861) q[1];
sx q[1];
rz(-1.4844196) q[1];
sx q[1];
rz(0.84815971) q[1];
rz(-pi) q[2];
rz(-0.72569287) q[3];
sx q[3];
rz(-1.4376493) q[3];
sx q[3];
rz(2.468688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1802804) q[2];
sx q[2];
rz(-1.140241) q[2];
sx q[2];
rz(0.22892924) q[2];
rz(1.3198352) q[3];
sx q[3];
rz(-1.6596551) q[3];
sx q[3];
rz(2.2975217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9127386) q[0];
sx q[0];
rz(-1.5874533) q[0];
sx q[0];
rz(-1.9629021) q[0];
rz(1.2294058) q[1];
sx q[1];
rz(-0.39437672) q[1];
sx q[1];
rz(1.2961402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79149198) q[0];
sx q[0];
rz(-1.7069567) q[0];
sx q[0];
rz(-1.8541992) q[0];
rz(-0.77028124) q[2];
sx q[2];
rz(-0.89777032) q[2];
sx q[2];
rz(1.9865303) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7414822) q[1];
sx q[1];
rz(-2.0026836) q[1];
sx q[1];
rz(-0.53439616) q[1];
x q[2];
rz(-1.6722176) q[3];
sx q[3];
rz(-1.0857255) q[3];
sx q[3];
rz(-0.17899638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.487315) q[2];
sx q[2];
rz(-0.9026022) q[2];
sx q[2];
rz(-2.8793092) q[2];
rz(-2.8413963) q[3];
sx q[3];
rz(-1.9191977) q[3];
sx q[3];
rz(0.20849553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7009785) q[0];
sx q[0];
rz(-2.746026) q[0];
sx q[0];
rz(1.3025008) q[0];
rz(2.473096) q[1];
sx q[1];
rz(-1.3553456) q[1];
sx q[1];
rz(2.1065333) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44159206) q[0];
sx q[0];
rz(-1.3574575) q[0];
sx q[0];
rz(1.160303) q[0];
rz(-pi) q[1];
rz(0.18685734) q[2];
sx q[2];
rz(-1.9344988) q[2];
sx q[2];
rz(-2.0314558) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0332843) q[1];
sx q[1];
rz(-2.7893745) q[1];
sx q[1];
rz(-1.3501106) q[1];
rz(-0.97888246) q[3];
sx q[3];
rz(-0.91875263) q[3];
sx q[3];
rz(-2.121821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11789007) q[2];
sx q[2];
rz(-2.1390476) q[2];
sx q[2];
rz(1.4855851) q[2];
rz(-2.8591136) q[3];
sx q[3];
rz(-1.3586724) q[3];
sx q[3];
rz(0.43454596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78541237) q[0];
sx q[0];
rz(-1.0419351) q[0];
sx q[0];
rz(2.8670512) q[0];
rz(-1.2046332) q[1];
sx q[1];
rz(-0.58012539) q[1];
sx q[1];
rz(0.56142941) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0208541) q[0];
sx q[0];
rz(-2.0971813) q[0];
sx q[0];
rz(0.83197439) q[0];
x q[1];
rz(-2.4653685) q[2];
sx q[2];
rz(-1.8100139) q[2];
sx q[2];
rz(-2.4382811) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9778487) q[1];
sx q[1];
rz(-1.1625966) q[1];
sx q[1];
rz(0.331649) q[1];
rz(-pi) q[2];
rz(0.26430126) q[3];
sx q[3];
rz(-1.9076548) q[3];
sx q[3];
rz(-1.5764396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5795035) q[2];
sx q[2];
rz(-2.0359437) q[2];
sx q[2];
rz(-1.4037464) q[2];
rz(1.3459407) q[3];
sx q[3];
rz(-2.3965049) q[3];
sx q[3];
rz(2.6534206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3812934) q[0];
sx q[0];
rz(-2.5385222) q[0];
sx q[0];
rz(2.2316933) q[0];
rz(-1.1970041) q[1];
sx q[1];
rz(-2.0884114) q[1];
sx q[1];
rz(-0.88814703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040374856) q[0];
sx q[0];
rz(-1.4618317) q[0];
sx q[0];
rz(3.0473112) q[0];
rz(-pi) q[1];
rz(-1.5106701) q[2];
sx q[2];
rz(-2.2297568) q[2];
sx q[2];
rz(1.2252285) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3067813) q[1];
sx q[1];
rz(-1.957292) q[1];
sx q[1];
rz(-0.67910925) q[1];
rz(-pi) q[2];
rz(3.1040583) q[3];
sx q[3];
rz(-1.8662631) q[3];
sx q[3];
rz(-0.42643828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1278648) q[2];
sx q[2];
rz(-1.4415386) q[2];
sx q[2];
rz(-2.0590651) q[2];
rz(-2.9108289) q[3];
sx q[3];
rz(-1.328238) q[3];
sx q[3];
rz(2.6575991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8312663) q[0];
sx q[0];
rz(-0.75208298) q[0];
sx q[0];
rz(0.58151522) q[0];
rz(-0.61559081) q[1];
sx q[1];
rz(-2.1938727) q[1];
sx q[1];
rz(1.8448578) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.588284) q[0];
sx q[0];
rz(-3.1206312) q[0];
sx q[0];
rz(-2.0901438) q[0];
x q[1];
rz(2.7974772) q[2];
sx q[2];
rz(-1.5798347) q[2];
sx q[2];
rz(-1.8923541) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.91528581) q[1];
sx q[1];
rz(-1.9939594) q[1];
sx q[1];
rz(-0.28171119) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1704135) q[3];
sx q[3];
rz(-0.93850157) q[3];
sx q[3];
rz(2.1454772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7221308) q[2];
sx q[2];
rz(-0.1408793) q[2];
sx q[2];
rz(0.44931832) q[2];
rz(2.8214473) q[3];
sx q[3];
rz(-1.3195427) q[3];
sx q[3];
rz(-1.8951353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.768059) q[0];
sx q[0];
rz(-2.4644485) q[0];
sx q[0];
rz(1.2039536) q[0];
rz(1.3846579) q[1];
sx q[1];
rz(-0.23778267) q[1];
sx q[1];
rz(2.3429088) q[1];
rz(2.8821386) q[2];
sx q[2];
rz(-0.96895192) q[2];
sx q[2];
rz(-0.12018991) q[2];
rz(2.0417236) q[3];
sx q[3];
rz(-1.2237329) q[3];
sx q[3];
rz(1.208186) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
