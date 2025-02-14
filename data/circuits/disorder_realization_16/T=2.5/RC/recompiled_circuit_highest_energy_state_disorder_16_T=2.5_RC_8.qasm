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
rz(0.056869153) q[0];
sx q[0];
rz(-0.19357227) q[0];
sx q[0];
rz(-0.40785664) q[0];
rz(3.1240533) q[1];
sx q[1];
rz(-1.2516302) q[1];
sx q[1];
rz(1.6012021) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33690573) q[0];
sx q[0];
rz(-2.0077472) q[0];
sx q[0];
rz(3.0880465) q[0];
rz(0.45405252) q[2];
sx q[2];
rz(-1.4489929) q[2];
sx q[2];
rz(1.4920838) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1931128) q[1];
sx q[1];
rz(-2.6958637) q[1];
sx q[1];
rz(2.6325001) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61397378) q[3];
sx q[3];
rz(-2.2212265) q[3];
sx q[3];
rz(2.2681595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5876329) q[2];
sx q[2];
rz(-1.0593869) q[2];
sx q[2];
rz(2.9239192) q[2];
rz(2.830128) q[3];
sx q[3];
rz(-2.5296827) q[3];
sx q[3];
rz(-1.1658839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72164732) q[0];
sx q[0];
rz(-2.8657275) q[0];
sx q[0];
rz(0.69197792) q[0];
rz(-2.3852589) q[1];
sx q[1];
rz(-2.6192009) q[1];
sx q[1];
rz(1.5914894) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78849906) q[0];
sx q[0];
rz(-2.278557) q[0];
sx q[0];
rz(-1.4097296) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1728376) q[2];
sx q[2];
rz(-1.6670456) q[2];
sx q[2];
rz(0.41589662) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8263232) q[1];
sx q[1];
rz(-0.20139748) q[1];
sx q[1];
rz(-2.4884495) q[1];
rz(1.0106929) q[3];
sx q[3];
rz(-1.1385131) q[3];
sx q[3];
rz(3.0280857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1456566) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(0.23925979) q[2];
rz(-1.4907106) q[3];
sx q[3];
rz(-1.2942634) q[3];
sx q[3];
rz(1.2283121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8656798) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(-0.014634125) q[0];
rz(0.89302653) q[1];
sx q[1];
rz(-2.1664186) q[1];
sx q[1];
rz(-0.21569529) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6071226) q[0];
sx q[0];
rz(-1.5661582) q[0];
sx q[0];
rz(-1.5195373) q[0];
rz(-0.65962445) q[2];
sx q[2];
rz(-1.8462034) q[2];
sx q[2];
rz(0.026247488) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0552935) q[1];
sx q[1];
rz(-2.8759416) q[1];
sx q[1];
rz(-2.6363027) q[1];
rz(2.2641597) q[3];
sx q[3];
rz(-1.0022396) q[3];
sx q[3];
rz(1.3662149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0821685) q[2];
sx q[2];
rz(-1.1801722) q[2];
sx q[2];
rz(-0.99349418) q[2];
rz(-0.70837402) q[3];
sx q[3];
rz(-1.1185442) q[3];
sx q[3];
rz(-0.12349252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3941536) q[0];
sx q[0];
rz(-0.49462947) q[0];
sx q[0];
rz(-2.4476449) q[0];
rz(2.8548062) q[1];
sx q[1];
rz(-1.6096121) q[1];
sx q[1];
rz(2.0538816) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3877849) q[0];
sx q[0];
rz(-1.0669363) q[0];
sx q[0];
rz(0.14150454) q[0];
rz(-pi) q[1];
rz(-1.0469646) q[2];
sx q[2];
rz(-2.0598542) q[2];
sx q[2];
rz(-0.2917052) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3301311) q[1];
sx q[1];
rz(-2.4566659) q[1];
sx q[1];
rz(-2.0622938) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4272653) q[3];
sx q[3];
rz(-1.5363524) q[3];
sx q[3];
rz(-1.2726651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22471681) q[2];
sx q[2];
rz(-1.1563053) q[2];
sx q[2];
rz(-1.3678331) q[2];
rz(1.2380838) q[3];
sx q[3];
rz(-1.5725458) q[3];
sx q[3];
rz(0.21627538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5424159) q[0];
sx q[0];
rz(-1.8734064) q[0];
sx q[0];
rz(0.61781484) q[0];
rz(2.0333596) q[1];
sx q[1];
rz(-2.8384659) q[1];
sx q[1];
rz(-0.88315001) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9734479) q[0];
sx q[0];
rz(-0.75605481) q[0];
sx q[0];
rz(0.010764695) q[0];
rz(0.3381392) q[2];
sx q[2];
rz(-2.6827742) q[2];
sx q[2];
rz(0.79566075) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2715192) q[1];
sx q[1];
rz(-0.49527676) q[1];
sx q[1];
rz(-1.7090767) q[1];
rz(-pi) q[2];
rz(2.1843076) q[3];
sx q[3];
rz(-0.64910474) q[3];
sx q[3];
rz(1.7136991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9250179) q[2];
sx q[2];
rz(-0.25187945) q[2];
sx q[2];
rz(1.3612755) q[2];
rz(1.8682293) q[3];
sx q[3];
rz(-1.7073771) q[3];
sx q[3];
rz(0.98289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044416044) q[0];
sx q[0];
rz(-1.8955078) q[0];
sx q[0];
rz(0.62028766) q[0];
rz(-0.39707956) q[1];
sx q[1];
rz(-0.47080165) q[1];
sx q[1];
rz(1.9420067) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7714149) q[0];
sx q[0];
rz(-1.4645029) q[0];
sx q[0];
rz(0.021076963) q[0];
x q[1];
rz(-2.1753005) q[2];
sx q[2];
rz(-2.8859865) q[2];
sx q[2];
rz(-2.5641455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5127677) q[1];
sx q[1];
rz(-0.56124748) q[1];
sx q[1];
rz(2.7732549) q[1];
rz(-pi) q[2];
rz(0.34440094) q[3];
sx q[3];
rz(-1.6823497) q[3];
sx q[3];
rz(-1.2408181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.24641307) q[2];
sx q[2];
rz(-1.2440888) q[2];
sx q[2];
rz(-2.4708774) q[2];
rz(-0.29117584) q[3];
sx q[3];
rz(-0.61101919) q[3];
sx q[3];
rz(-2.5197855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7156242) q[0];
sx q[0];
rz(-1.7142897) q[0];
sx q[0];
rz(-0.17844644) q[0];
rz(0.87002358) q[1];
sx q[1];
rz(-0.88302892) q[1];
sx q[1];
rz(-2.3387486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4492396) q[0];
sx q[0];
rz(-0.62004161) q[0];
sx q[0];
rz(-3.0983144) q[0];
x q[1];
rz(2.5576791) q[2];
sx q[2];
rz(-1.7758992) q[2];
sx q[2];
rz(0.55196229) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72910129) q[1];
sx q[1];
rz(-1.2454528) q[1];
sx q[1];
rz(0.71112432) q[1];
rz(-pi) q[2];
rz(-2.1986897) q[3];
sx q[3];
rz(-1.781309) q[3];
sx q[3];
rz(2.877416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7650083) q[2];
sx q[2];
rz(-1.2976126) q[2];
sx q[2];
rz(0.89361781) q[2];
rz(1.5417967) q[3];
sx q[3];
rz(-1.6518281) q[3];
sx q[3];
rz(0.60683513) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1641418) q[0];
sx q[0];
rz(-2.7315388) q[0];
sx q[0];
rz(-2.9264911) q[0];
rz(-2.2970301) q[1];
sx q[1];
rz(-1.3203878) q[1];
sx q[1];
rz(-2.3473158) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2820123) q[0];
sx q[0];
rz(-1.5589801) q[0];
sx q[0];
rz(-0.026547308) q[0];
rz(-2.568001) q[2];
sx q[2];
rz(-0.52917889) q[2];
sx q[2];
rz(2.0647788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.46103288) q[1];
sx q[1];
rz(-2.5941412) q[1];
sx q[1];
rz(1.3957109) q[1];
x q[2];
rz(-1.896554) q[3];
sx q[3];
rz(-2.3255201) q[3];
sx q[3];
rz(-1.9297502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8636785) q[2];
sx q[2];
rz(-1.9387551) q[2];
sx q[2];
rz(1.7318783) q[2];
rz(0.75285161) q[3];
sx q[3];
rz(-1.1466305) q[3];
sx q[3];
rz(-2.1394155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3442605) q[0];
sx q[0];
rz(-0.40818885) q[0];
sx q[0];
rz(-2.1687188) q[0];
rz(1.5219888) q[1];
sx q[1];
rz(-1.3643967) q[1];
sx q[1];
rz(-2.5111228) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79463835) q[0];
sx q[0];
rz(-1.7525867) q[0];
sx q[0];
rz(1.4304881) q[0];
rz(-pi) q[1];
rz(-1.5233598) q[2];
sx q[2];
rz(-2.9844797) q[2];
sx q[2];
rz(-2.4542798) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9520397) q[1];
sx q[1];
rz(-1.5209043) q[1];
sx q[1];
rz(-2.3734809) q[1];
rz(-0.34548605) q[3];
sx q[3];
rz(-2.957323) q[3];
sx q[3];
rz(2.1994182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1260599) q[2];
sx q[2];
rz(-1.9893179) q[2];
sx q[2];
rz(-1.3014334) q[2];
rz(-1.556501) q[3];
sx q[3];
rz(-2.3157178) q[3];
sx q[3];
rz(-0.41527709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5179317) q[0];
sx q[0];
rz(-2.9908337) q[0];
sx q[0];
rz(-2.6483722) q[0];
rz(-1.8863691) q[1];
sx q[1];
rz(-1.247765) q[1];
sx q[1];
rz(0.36466041) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4987558) q[0];
sx q[0];
rz(-0.36223534) q[0];
sx q[0];
rz(-1.1722159) q[0];
x q[1];
rz(-1.0027867) q[2];
sx q[2];
rz(-0.73633654) q[2];
sx q[2];
rz(1.9970837) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9182049) q[1];
sx q[1];
rz(-0.80092309) q[1];
sx q[1];
rz(-0.078322874) q[1];
x q[2];
rz(-3.10121) q[3];
sx q[3];
rz(-1.8869487) q[3];
sx q[3];
rz(-0.97163661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6898592) q[2];
sx q[2];
rz(-2.1351337) q[2];
sx q[2];
rz(-2.149392) q[2];
rz(-0.36710468) q[3];
sx q[3];
rz(-1.4398984) q[3];
sx q[3];
rz(-2.4632857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.622396) q[0];
sx q[0];
rz(-1.665103) q[0];
sx q[0];
rz(-1.7892224) q[0];
rz(2.2182111) q[1];
sx q[1];
rz(-2.2877749) q[1];
sx q[1];
rz(1.9705082) q[1];
rz(2.1869356) q[2];
sx q[2];
rz(-2.8980394) q[2];
sx q[2];
rz(0.63799636) q[2];
rz(-0.70960511) q[3];
sx q[3];
rz(-0.81021877) q[3];
sx q[3];
rz(2.6859663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
