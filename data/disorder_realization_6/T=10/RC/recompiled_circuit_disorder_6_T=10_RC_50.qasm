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
rz(-1.7237741) q[0];
sx q[0];
rz(-0.56086993) q[0];
rz(1.1129192) q[1];
sx q[1];
rz(-1.7634044) q[1];
sx q[1];
rz(1.2150432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5117447) q[0];
sx q[0];
rz(-1.4667604) q[0];
sx q[0];
rz(1.758979) q[0];
rz(-2.4807793) q[2];
sx q[2];
rz(-1.8006969) q[2];
sx q[2];
rz(1.4290609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.558555) q[1];
sx q[1];
rz(-1.8999294) q[1];
sx q[1];
rz(-2.1117044) q[1];
rz(1.3052985) q[3];
sx q[3];
rz(-1.4269097) q[3];
sx q[3];
rz(3.0778411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3540196) q[2];
sx q[2];
rz(-2.188787) q[2];
sx q[2];
rz(-0.18307486) q[2];
rz(2.7637774) q[3];
sx q[3];
rz(-2.0928045) q[3];
sx q[3];
rz(2.8474076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8437682) q[0];
sx q[0];
rz(-2.4968708) q[0];
sx q[0];
rz(0.077117292) q[0];
rz(-0.33879694) q[1];
sx q[1];
rz(-2.0270551) q[1];
sx q[1];
rz(-1.6024626) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3829271) q[0];
sx q[0];
rz(-0.59249741) q[0];
sx q[0];
rz(1.4325607) q[0];
rz(-2.9039731) q[2];
sx q[2];
rz(-2.3887206) q[2];
sx q[2];
rz(1.876229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96670818) q[1];
sx q[1];
rz(-2.0497353) q[1];
sx q[1];
rz(0.19660463) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7178571) q[3];
sx q[3];
rz(-1.7859922) q[3];
sx q[3];
rz(-2.7620897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2960647) q[2];
sx q[2];
rz(-1.277593) q[2];
sx q[2];
rz(-2.4831333) q[2];
rz(-2.9902839) q[3];
sx q[3];
rz(-2.1189809) q[3];
sx q[3];
rz(-2.4466799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028458683) q[0];
sx q[0];
rz(-2.3582393) q[0];
sx q[0];
rz(-0.43310305) q[0];
rz(1.1921047) q[1];
sx q[1];
rz(-1.2116218) q[1];
sx q[1];
rz(-0.55535299) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99291891) q[0];
sx q[0];
rz(-0.86704555) q[0];
sx q[0];
rz(-0.91627319) q[0];
rz(-2.0434415) q[2];
sx q[2];
rz(-0.61908365) q[2];
sx q[2];
rz(-1.667779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8086116) q[1];
sx q[1];
rz(-0.65578991) q[1];
sx q[1];
rz(-1.8379184) q[1];
x q[2];
rz(1.8954574) q[3];
sx q[3];
rz(-0.56926308) q[3];
sx q[3];
rz(3.0666921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73734036) q[2];
sx q[2];
rz(-0.78812391) q[2];
sx q[2];
rz(1.2505442) q[2];
rz(0.2441497) q[3];
sx q[3];
rz(-1.8593676) q[3];
sx q[3];
rz(1.6916493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26043949) q[0];
sx q[0];
rz(-0.45757159) q[0];
sx q[0];
rz(2.326791) q[0];
rz(-1.3793777) q[1];
sx q[1];
rz(-0.35019362) q[1];
sx q[1];
rz(0.25517685) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4677306) q[0];
sx q[0];
rz(-0.47463372) q[0];
sx q[0];
rz(-0.69068308) q[0];
x q[1];
rz(-2.5023979) q[2];
sx q[2];
rz(-1.8393469) q[2];
sx q[2];
rz(-0.07721363) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0604917) q[1];
sx q[1];
rz(-1.2730518) q[1];
sx q[1];
rz(-0.22025073) q[1];
rz(-0.051024036) q[3];
sx q[3];
rz(-1.050356) q[3];
sx q[3];
rz(-2.737962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2531551) q[2];
sx q[2];
rz(-1.5506813) q[2];
sx q[2];
rz(-2.9684084) q[2];
rz(-2.611768) q[3];
sx q[3];
rz(-0.14557043) q[3];
sx q[3];
rz(-0.1023275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2816876) q[0];
sx q[0];
rz(-1.4929993) q[0];
sx q[0];
rz(1.3758855) q[0];
rz(-1.2777404) q[1];
sx q[1];
rz(-0.81218305) q[1];
sx q[1];
rz(-3.0854991) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8431906) q[0];
sx q[0];
rz(-2.139233) q[0];
sx q[0];
rz(-2.0060904) q[0];
rz(-2.202583) q[2];
sx q[2];
rz(-0.54982215) q[2];
sx q[2];
rz(0.75013559) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6009112) q[1];
sx q[1];
rz(-1.8735421) q[1];
sx q[1];
rz(1.5860228) q[1];
x q[2];
rz(1.0104996) q[3];
sx q[3];
rz(-1.4951402) q[3];
sx q[3];
rz(1.2922985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.00099480199) q[2];
sx q[2];
rz(-2.2237015) q[2];
sx q[2];
rz(2.7094005) q[2];
rz(2.2473992) q[3];
sx q[3];
rz(-2.0420572) q[3];
sx q[3];
rz(-1.6754707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.4636661) q[1];
sx q[1];
rz(0.9544968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.944825) q[0];
sx q[0];
rz(-1.4155354) q[0];
sx q[0];
rz(-1.0374271) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6184094) q[2];
sx q[2];
rz(-2.5104668) q[2];
sx q[2];
rz(-2.7472251) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2973605) q[1];
sx q[1];
rz(-2.4095222) q[1];
sx q[1];
rz(-0.75433235) q[1];
x q[2];
rz(-2.9783863) q[3];
sx q[3];
rz(-1.8582134) q[3];
sx q[3];
rz(-0.67374574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.59297562) q[2];
sx q[2];
rz(-1.9081215) q[2];
sx q[2];
rz(-1.0423638) q[2];
rz(-2.7029165) q[3];
sx q[3];
rz(-1.0500267) q[3];
sx q[3];
rz(1.3180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577268) q[0];
sx q[0];
rz(-0.23290578) q[0];
sx q[0];
rz(0.74321157) q[0];
rz(1.5076393) q[1];
sx q[1];
rz(-2.4217024) q[1];
sx q[1];
rz(-0.61002237) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9305206) q[0];
sx q[0];
rz(-1.0809582) q[0];
sx q[0];
rz(1.8563489) q[0];
x q[1];
rz(2.2248613) q[2];
sx q[2];
rz(-1.6628633) q[2];
sx q[2];
rz(2.9031861) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94177946) q[1];
sx q[1];
rz(-2.4738414) q[1];
sx q[1];
rz(-1.6112531) q[1];
rz(-pi) q[2];
rz(1.0878629) q[3];
sx q[3];
rz(-0.42290877) q[3];
sx q[3];
rz(-1.203323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8053749) q[2];
sx q[2];
rz(-1.4423794) q[2];
sx q[2];
rz(-2.2231893) q[2];
rz(-1.5504799) q[3];
sx q[3];
rz(-0.95033002) q[3];
sx q[3];
rz(-2.7526855) q[3];
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
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.780705) q[0];
sx q[0];
rz(-2.4724859) q[0];
sx q[0];
rz(-1.6280744) q[0];
rz(-0.52945119) q[1];
sx q[1];
rz(-2.0748731) q[1];
sx q[1];
rz(2.4050074) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37305957) q[0];
sx q[0];
rz(-0.031950843) q[0];
sx q[0];
rz(1.2241227) q[0];
x q[1];
rz(3.0326764) q[2];
sx q[2];
rz(-1.250259) q[2];
sx q[2];
rz(-0.045217302) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.089162) q[1];
sx q[1];
rz(-0.4757291) q[1];
sx q[1];
rz(-2.9176941) q[1];
rz(-pi) q[2];
rz(2.1804817) q[3];
sx q[3];
rz(-2.612252) q[3];
sx q[3];
rz(1.5101658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1071876) q[2];
sx q[2];
rz(-1.2074869) q[2];
sx q[2];
rz(-2.4576808) q[2];
rz(1.9125787) q[3];
sx q[3];
rz(-1.3701655) q[3];
sx q[3];
rz(1.3945403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443611) q[0];
sx q[0];
rz(-1.5988388) q[0];
sx q[0];
rz(0.18572447) q[0];
rz(-2.1445403) q[1];
sx q[1];
rz(-1.2652218) q[1];
sx q[1];
rz(-0.7448147) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1574402) q[0];
sx q[0];
rz(-0.8343578) q[0];
sx q[0];
rz(-2.0122583) q[0];
x q[1];
rz(1.1698193) q[2];
sx q[2];
rz(-2.3224761) q[2];
sx q[2];
rz(-0.70105201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.71036584) q[1];
sx q[1];
rz(-1.2214298) q[1];
sx q[1];
rz(2.9037895) q[1];
x q[2];
rz(-1.1780147) q[3];
sx q[3];
rz(-0.82735705) q[3];
sx q[3];
rz(-2.3208997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0361438) q[2];
sx q[2];
rz(-0.75573409) q[2];
sx q[2];
rz(-1.194681) q[2];
rz(-2.1448994) q[3];
sx q[3];
rz(-1.2160622) q[3];
sx q[3];
rz(0.99635807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74334082) q[0];
sx q[0];
rz(-2.353459) q[0];
sx q[0];
rz(-0.40400305) q[0];
rz(0.031127302) q[1];
sx q[1];
rz(-1.6571836) q[1];
sx q[1];
rz(-1.1709447) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0572646) q[0];
sx q[0];
rz(-1.9976166) q[0];
sx q[0];
rz(1.5469993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77565907) q[2];
sx q[2];
rz(-1.7663029) q[2];
sx q[2];
rz(1.0722216) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7217692) q[1];
sx q[1];
rz(-3.1173919) q[1];
sx q[1];
rz(1.9355378) q[1];
rz(-pi) q[2];
rz(0.42697866) q[3];
sx q[3];
rz(-1.3037762) q[3];
sx q[3];
rz(-0.094735183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6955473) q[2];
sx q[2];
rz(-1.3617159) q[2];
sx q[2];
rz(-2.5496303) q[2];
rz(-0.56636089) q[3];
sx q[3];
rz(-2.9768894) q[3];
sx q[3];
rz(1.5238354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3175209) q[0];
sx q[0];
rz(-2.1614647) q[0];
sx q[0];
rz(1.9807057) q[0];
rz(3.042165) q[1];
sx q[1];
rz(-1.2482523) q[1];
sx q[1];
rz(-2.0773239) q[1];
rz(-0.912491) q[2];
sx q[2];
rz(-2.1952663) q[2];
sx q[2];
rz(-1.2637539) q[2];
rz(-3.1254461) q[3];
sx q[3];
rz(-1.9042249) q[3];
sx q[3];
rz(-0.92845542) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];