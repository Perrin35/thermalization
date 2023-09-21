OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.52580994) q[0];
sx q[0];
rz(4.5594112) q[0];
sx q[0];
rz(8.863908) q[0];
rz(-2.0286735) q[1];
sx q[1];
rz(-1.3781883) q[1];
sx q[1];
rz(1.9265494) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44170609) q[0];
sx q[0];
rz(-0.21472782) q[0];
sx q[0];
rz(1.0617274) q[0];
rz(-pi) q[1];
rz(-2.7772929) q[2];
sx q[2];
rz(-0.69395739) q[2];
sx q[2];
rz(-2.7147164) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3455968) q[1];
sx q[1];
rz(-2.0797634) q[1];
sx q[1];
rz(-2.7624346) q[1];
x q[2];
rz(-2.9925572) q[3];
sx q[3];
rz(-1.8334853) q[3];
sx q[3];
rz(-1.5460154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.78757301) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(0.18307486) q[2];
rz(0.37781528) q[3];
sx q[3];
rz(-1.0487882) q[3];
sx q[3];
rz(2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29782444) q[0];
sx q[0];
rz(-0.64472187) q[0];
sx q[0];
rz(-0.077117292) q[0];
rz(0.33879694) q[1];
sx q[1];
rz(-1.1145376) q[1];
sx q[1];
rz(1.5391301) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7586655) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(-1.4325607) q[0];
rz(-1.7878754) q[2];
sx q[2];
rz(-2.2976544) q[2];
sx q[2];
rz(2.196687) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.96670818) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(0.19660463) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3354524) q[3];
sx q[3];
rz(-1.1574405) q[3];
sx q[3];
rz(1.2873161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2960647) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-0.65845931) q[2];
rz(-0.15130875) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.113134) q[0];
sx q[0];
rz(-0.78335339) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.9299709) q[1];
sx q[1];
rz(0.55535299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0387602) q[0];
sx q[0];
rz(-1.0881249) q[0];
sx q[0];
rz(-2.32248) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31366445) q[2];
sx q[2];
rz(-1.0278388) q[2];
sx q[2];
rz(-1.107159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8086116) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(-1.8379184) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2461353) q[3];
sx q[3];
rz(-2.5723296) q[3];
sx q[3];
rz(-0.07490052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(-1.8910485) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(-1.4499433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(-2.326791) q[0];
rz(1.762215) q[1];
sx q[1];
rz(-2.791399) q[1];
sx q[1];
rz(-0.25517685) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92514738) q[0];
sx q[0];
rz(-1.9307923) q[0];
sx q[0];
rz(1.2544592) q[0];
rz(-0.63919477) q[2];
sx q[2];
rz(-1.3022458) q[2];
sx q[2];
rz(-0.07721363) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7123588) q[1];
sx q[1];
rz(-0.36839596) q[1];
sx q[1];
rz(2.1894987) q[1];
rz(-pi) q[2];
rz(3.0905686) q[3];
sx q[3];
rz(-2.0912366) q[3];
sx q[3];
rz(-0.40363064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8884376) q[2];
sx q[2];
rz(-1.5909114) q[2];
sx q[2];
rz(0.17318428) q[2];
rz(-0.52982461) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.6485933) q[0];
sx q[0];
rz(1.7657071) q[0];
rz(-1.2777404) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(-3.0854991) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027095196) q[0];
sx q[0];
rz(-1.9341015) q[0];
sx q[0];
rz(-0.61371213) q[0];
rz(-pi) q[1];
rz(0.9390097) q[2];
sx q[2];
rz(-2.5917705) q[2];
sx q[2];
rz(-0.75013559) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6009112) q[1];
sx q[1];
rz(-1.2680506) q[1];
sx q[1];
rz(-1.5555698) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4291184) q[3];
sx q[3];
rz(-2.5767527) q[3];
sx q[3];
rz(2.7431938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.00099480199) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(-2.7094005) q[2];
rz(2.2473992) q[3];
sx q[3];
rz(-1.0995355) q[3];
sx q[3];
rz(1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0284001) q[0];
sx q[0];
rz(-2.2560461) q[0];
sx q[0];
rz(0.64754852) q[0];
rz(-1.8796857) q[1];
sx q[1];
rz(-1.6779265) q[1];
sx q[1];
rz(-0.9544968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19676767) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(1.0374271) q[0];
rz(-pi) q[1];
rz(-3.1068222) q[2];
sx q[2];
rz(-2.2010942) q[2];
sx q[2];
rz(2.6882753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1986188) q[1];
sx q[1];
rz(-2.0795515) q[1];
sx q[1];
rz(1.0191304) q[1];
rz(-pi) q[2];
rz(1.2797221) q[3];
sx q[3];
rz(-1.4143412) q[3];
sx q[3];
rz(2.1978956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.548617) q[2];
sx q[2];
rz(-1.2334712) q[2];
sx q[2];
rz(-1.0423638) q[2];
rz(-0.43867612) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.8235122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(2.3983811) q[0];
rz(1.6339533) q[1];
sx q[1];
rz(-0.71989027) q[1];
sx q[1];
rz(2.5315703) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7689432) q[0];
sx q[0];
rz(-2.5805051) q[0];
sx q[0];
rz(2.6555496) q[0];
x q[1];
rz(0.91673135) q[2];
sx q[2];
rz(-1.4787294) q[2];
sx q[2];
rz(-0.23840657) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.94177946) q[1];
sx q[1];
rz(-2.4738414) q[1];
sx q[1];
rz(-1.6112531) q[1];
x q[2];
rz(2.9355572) q[3];
sx q[3];
rz(-1.1987975) q[3];
sx q[3];
rz(0.68148617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.33621776) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(-1.5911128) q[3];
sx q[3];
rz(-2.1912626) q[3];
sx q[3];
rz(-2.7526855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36088762) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(2.6121415) q[1];
sx q[1];
rz(-1.0667195) q[1];
sx q[1];
rz(-2.4050074) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37305957) q[0];
sx q[0];
rz(-3.1096418) q[0];
sx q[0];
rz(1.2241227) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2484776) q[2];
sx q[2];
rz(-1.6741447) q[2];
sx q[2];
rz(-1.5600187) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3031993) q[1];
sx q[1];
rz(-2.0337078) q[1];
sx q[1];
rz(-1.4569015) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96111091) q[3];
sx q[3];
rz(-0.52934066) q[3];
sx q[3];
rz(1.5101658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0344051) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(-2.4576808) q[2];
rz(-1.2290139) q[3];
sx q[3];
rz(-1.7714272) q[3];
sx q[3];
rz(1.7470523) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1972315) q[0];
sx q[0];
rz(-1.5427538) q[0];
sx q[0];
rz(-0.18572447) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.8763708) q[1];
sx q[1];
rz(-2.396778) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72070044) q[0];
sx q[0];
rz(-1.2487131) q[0];
sx q[0];
rz(-2.3548404) q[0];
rz(-1.1698193) q[2];
sx q[2];
rz(-0.81911659) q[2];
sx q[2];
rz(-0.70105201) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71036584) q[1];
sx q[1];
rz(-1.9201628) q[1];
sx q[1];
rz(-2.9037895) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9635779) q[3];
sx q[3];
rz(-0.82735705) q[3];
sx q[3];
rz(0.82069293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-2.3858586) q[2];
sx q[2];
rz(1.194681) q[2];
rz(0.99669325) q[3];
sx q[3];
rz(-1.9255305) q[3];
sx q[3];
rz(-0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-2.3982518) q[0];
sx q[0];
rz(-0.78813362) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(-3.1104654) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.1709447) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4963213) q[0];
sx q[0];
rz(-1.5924581) q[0];
sx q[0];
rz(-2.7146656) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8413713) q[2];
sx q[2];
rz(-2.3279394) q[2];
sx q[2];
rz(-2.4547581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.4198235) q[1];
sx q[1];
rz(-0.024200736) q[1];
sx q[1];
rz(-1.2060549) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42697866) q[3];
sx q[3];
rz(-1.8378165) q[3];
sx q[3];
rz(-3.0468575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4460454) q[2];
sx q[2];
rz(-1.7798767) q[2];
sx q[2];
rz(-2.5496303) q[2];
rz(-2.5752318) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(-1.5238354) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3175209) q[0];
sx q[0];
rz(-0.98012797) q[0];
sx q[0];
rz(-1.160887) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(2.4026985) q[2];
sx q[2];
rz(-2.0901491) q[2];
sx q[2];
rz(0.73170589) q[2];
rz(1.6173784) q[3];
sx q[3];
rz(-0.33380476) q[3];
sx q[3];
rz(-0.87915626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];