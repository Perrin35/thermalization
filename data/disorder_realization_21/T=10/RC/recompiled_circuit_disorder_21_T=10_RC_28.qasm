OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0857467) q[0];
sx q[0];
rz(-0.081781713) q[0];
sx q[0];
rz(-2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6591588) q[0];
sx q[0];
rz(-1.6560439) q[0];
sx q[0];
rz(1.2486588) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6104923) q[2];
sx q[2];
rz(-0.98449003) q[2];
sx q[2];
rz(2.5551978) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1134125) q[1];
sx q[1];
rz(-2.1537188) q[1];
sx q[1];
rz(-2.7229573) q[1];
rz(-pi) q[2];
rz(2.206771) q[3];
sx q[3];
rz(-1.9536195) q[3];
sx q[3];
rz(-2.8179907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6364608) q[2];
sx q[2];
rz(-0.59288609) q[2];
sx q[2];
rz(0.55603975) q[2];
rz(0.83267823) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(-2.1957943) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(-0.15727501) q[0];
rz(-2.8804624) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(-0.10903407) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1125688) q[0];
sx q[0];
rz(-1.3431664) q[0];
sx q[0];
rz(-0.17507041) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7686339) q[2];
sx q[2];
rz(-1.3126144) q[2];
sx q[2];
rz(2.3193662) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2687159) q[1];
sx q[1];
rz(-1.1040338) q[1];
sx q[1];
rz(-0.66653911) q[1];
rz(2.6208932) q[3];
sx q[3];
rz(-1.4200913) q[3];
sx q[3];
rz(0.13651383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9033501) q[2];
sx q[2];
rz(-1.1652596) q[2];
sx q[2];
rz(1.2634574) q[2];
rz(-0.3271099) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(-1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064421244) q[0];
sx q[0];
rz(-0.049296878) q[0];
sx q[0];
rz(-1.7984614) q[0];
rz(-2.893977) q[1];
sx q[1];
rz(-2.394948) q[1];
sx q[1];
rz(-2.6599191) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4108698) q[0];
sx q[0];
rz(-2.0559089) q[0];
sx q[0];
rz(-2.7184125) q[0];
rz(0.33294296) q[2];
sx q[2];
rz(-0.98368401) q[2];
sx q[2];
rz(-1.743403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18322769) q[1];
sx q[1];
rz(-1.0291568) q[1];
sx q[1];
rz(0.72802131) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89110156) q[3];
sx q[3];
rz(-2.292233) q[3];
sx q[3];
rz(-2.6654055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8032916) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(-2.6417007) q[2];
rz(0.56097427) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(-1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9445779) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(-2.5894077) q[0];
rz(1.588297) q[1];
sx q[1];
rz(-2.242656) q[1];
sx q[1];
rz(-1.8968556) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3361928) q[0];
sx q[0];
rz(-1.5691783) q[0];
sx q[0];
rz(-2.8021115) q[0];
rz(1.3696026) q[2];
sx q[2];
rz(-1.4428713) q[2];
sx q[2];
rz(0.97845562) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4501805) q[1];
sx q[1];
rz(-2.1069063) q[1];
sx q[1];
rz(0.24713534) q[1];
rz(2.5707385) q[3];
sx q[3];
rz(-1.7613162) q[3];
sx q[3];
rz(1.1754456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(-1.1506895) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(-2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078159049) q[0];
sx q[0];
rz(-2.3796191) q[0];
sx q[0];
rz(3.0601236) q[0];
rz(-0.062462656) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(-1.5030456) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.016184729) q[0];
sx q[0];
rz(-0.59690079) q[0];
sx q[0];
rz(1.3422658) q[0];
rz(-pi) q[1];
x q[1];
rz(2.224515) q[2];
sx q[2];
rz(-0.13609016) q[2];
sx q[2];
rz(-1.2166785) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11322184) q[1];
sx q[1];
rz(-1.8824667) q[1];
sx q[1];
rz(3.034364) q[1];
rz(-pi) q[2];
rz(1.2522069) q[3];
sx q[3];
rz(-1.0831523) q[3];
sx q[3];
rz(-2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(-1.2379237) q[2];
rz(-1.1226908) q[3];
sx q[3];
rz(-0.676238) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5313107) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(0.26671985) q[0];
rz(2.5807014) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(-0.7985324) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079433867) q[0];
sx q[0];
rz(-1.2957797) q[0];
sx q[0];
rz(-2.9955203) q[0];
rz(-pi) q[1];
rz(-1.2217667) q[2];
sx q[2];
rz(-1.0997699) q[2];
sx q[2];
rz(0.40086056) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16520271) q[1];
sx q[1];
rz(-0.573728) q[1];
sx q[1];
rz(-1.6673253) q[1];
rz(-2.5391606) q[3];
sx q[3];
rz(-2.2450387) q[3];
sx q[3];
rz(0.34365052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(0.36671656) q[2];
rz(1.2612873) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(-3.0453851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40903184) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(-0.60638705) q[0];
rz(-0.19730332) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(2.6775449) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1394135) q[0];
sx q[0];
rz(-2.3521949) q[0];
sx q[0];
rz(-2.1887357) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1214141) q[2];
sx q[2];
rz(-1.3300465) q[2];
sx q[2];
rz(-1.9764331) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5210919) q[1];
sx q[1];
rz(-1.8720086) q[1];
sx q[1];
rz(-1.0628994) q[1];
rz(-pi) q[2];
rz(-2.9052827) q[3];
sx q[3];
rz(-1.4910306) q[3];
sx q[3];
rz(0.41111708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2074034) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(-0.25804538) q[2];
rz(1.1856273) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946452) q[0];
sx q[0];
rz(-1.2807245) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(-3.0463468) q[1];
sx q[1];
rz(-2.1691599) q[1];
sx q[1];
rz(-1.4415178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7521116) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(1.585929) q[0];
x q[1];
rz(3.0870073) q[2];
sx q[2];
rz(-0.61331257) q[2];
sx q[2];
rz(-1.7577946) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6724832) q[1];
sx q[1];
rz(-1.4454578) q[1];
sx q[1];
rz(-3.010716) q[1];
x q[2];
rz(0.33393319) q[3];
sx q[3];
rz(-1.5302368) q[3];
sx q[3];
rz(3.0143152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9986481) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(-0.22658919) q[2];
rz(-2.6930124) q[3];
sx q[3];
rz(-1.6058763) q[3];
sx q[3];
rz(-2.3118238) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3806234) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(-1.3990078) q[0];
rz(0.31708583) q[1];
sx q[1];
rz(-1.4750907) q[1];
sx q[1];
rz(2.1549966) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6331659) q[0];
sx q[0];
rz(-1.0796483) q[0];
sx q[0];
rz(-1.3915865) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6096452) q[2];
sx q[2];
rz(-0.94594687) q[2];
sx q[2];
rz(0.22235409) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66227312) q[1];
sx q[1];
rz(-2.5179177) q[1];
sx q[1];
rz(0.24346607) q[1];
rz(-pi) q[2];
rz(-1.9429728) q[3];
sx q[3];
rz(-2.6937006) q[3];
sx q[3];
rz(0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(-0.26327291) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(2.6221361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(-1.7364527) q[0];
rz(-0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(-1.4155037) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6575359) q[0];
sx q[0];
rz(-0.41175479) q[0];
sx q[0];
rz(2.5176237) q[0];
rz(-pi) q[1];
rz(-2.892898) q[2];
sx q[2];
rz(-1.9242312) q[2];
sx q[2];
rz(0.24662185) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.532383) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(-1.3551559) q[1];
rz(-pi) q[2];
rz(-1.2763001) q[3];
sx q[3];
rz(-2.6101972) q[3];
sx q[3];
rz(-2.6607799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(0.12410513) q[2];
rz(2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(2.4805099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.538095) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(-0.30766906) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(-1.7514501) q[2];
sx q[2];
rz(-1.3144819) q[2];
sx q[2];
rz(-1.1983295) q[2];
rz(2.5849871) q[3];
sx q[3];
rz(-1.442933) q[3];
sx q[3];
rz(1.8419151) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
