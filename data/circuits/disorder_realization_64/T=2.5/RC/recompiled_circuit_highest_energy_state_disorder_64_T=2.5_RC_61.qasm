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
rz(-0.023802726) q[0];
sx q[0];
rz(4.2476141) q[0];
sx q[0];
rz(10.201693) q[0];
rz(-1.7493526) q[1];
sx q[1];
rz(-1.8267781) q[1];
sx q[1];
rz(-2.1652752) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38589729) q[0];
sx q[0];
rz(-1.160826) q[0];
sx q[0];
rz(-0.14571054) q[0];
rz(0.56371477) q[2];
sx q[2];
rz(-1.6275121) q[2];
sx q[2];
rz(0.67981718) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0529671) q[1];
sx q[1];
rz(-1.8421116) q[1];
sx q[1];
rz(-2.0120828) q[1];
rz(-0.99774811) q[3];
sx q[3];
rz(-1.5437417) q[3];
sx q[3];
rz(1.7535576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6355847) q[2];
sx q[2];
rz(-3.0878461) q[2];
sx q[2];
rz(-0.40679833) q[2];
rz(2.9721416) q[3];
sx q[3];
rz(-2.6120766) q[3];
sx q[3];
rz(1.0725526) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605543) q[0];
sx q[0];
rz(-2.9091703) q[0];
sx q[0];
rz(0.0037923092) q[0];
rz(0.077839851) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(-2.8357764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7821976) q[0];
sx q[0];
rz(-2.5919624) q[0];
sx q[0];
rz(-1.2122985) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3897116) q[2];
sx q[2];
rz(-1.4724178) q[2];
sx q[2];
rz(1.6461314) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0111573) q[1];
sx q[1];
rz(-0.88825916) q[1];
sx q[1];
rz(2.621317) q[1];
x q[2];
rz(-0.060896994) q[3];
sx q[3];
rz(-1.6346667) q[3];
sx q[3];
rz(-2.8012432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.99016142) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(2.4988417) q[2];
rz(-0.07240545) q[3];
sx q[3];
rz(-2.0824771) q[3];
sx q[3];
rz(-1.7806627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0533503) q[0];
sx q[0];
rz(-2.998816) q[0];
sx q[0];
rz(2.9852168) q[0];
rz(-0.0414255) q[1];
sx q[1];
rz(-2.5138469) q[1];
sx q[1];
rz(-1.590439) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56083) q[0];
sx q[0];
rz(-1.7770045) q[0];
sx q[0];
rz(-2.5872562) q[0];
rz(-1.5029491) q[2];
sx q[2];
rz(-2.2257651) q[2];
sx q[2];
rz(2.613435) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9030021) q[1];
sx q[1];
rz(-2.6005473) q[1];
sx q[1];
rz(2.9220086) q[1];
rz(0.66371347) q[3];
sx q[3];
rz(-0.98582375) q[3];
sx q[3];
rz(-3.0970517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40858832) q[2];
sx q[2];
rz(-2.6027347) q[2];
sx q[2];
rz(3.1206701) q[2];
rz(2.9544592) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(-3.0278897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1784096) q[0];
sx q[0];
rz(-2.8091176) q[0];
sx q[0];
rz(0.43854976) q[0];
rz(-1.5248388) q[1];
sx q[1];
rz(-2.8068145) q[1];
sx q[1];
rz(2.8964892) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49798508) q[0];
sx q[0];
rz(-1.4873355) q[0];
sx q[0];
rz(0.85026922) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7510277) q[2];
sx q[2];
rz(-0.3949983) q[2];
sx q[2];
rz(-2.4633138) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.72232258) q[1];
sx q[1];
rz(-0.55907226) q[1];
sx q[1];
rz(-1.30617) q[1];
rz(-1.2538337) q[3];
sx q[3];
rz(-0.84268236) q[3];
sx q[3];
rz(1.1945981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(-0.33176804) q[2];
rz(0.48745421) q[3];
sx q[3];
rz(-1.0468227) q[3];
sx q[3];
rz(0.99307466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29193923) q[0];
sx q[0];
rz(-1.4537469) q[0];
sx q[0];
rz(-2.3680903) q[0];
rz(-1.1812814) q[1];
sx q[1];
rz(-2.9995194) q[1];
sx q[1];
rz(-1.389651) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31598241) q[0];
sx q[0];
rz(-1.1706691) q[0];
sx q[0];
rz(1.1745625) q[0];
x q[1];
rz(-0.59771363) q[2];
sx q[2];
rz(-2.1897912) q[2];
sx q[2];
rz(-2.5368382) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4378499) q[1];
sx q[1];
rz(-1.653076) q[1];
sx q[1];
rz(2.6668496) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64719154) q[3];
sx q[3];
rz(-0.29509896) q[3];
sx q[3];
rz(1.5403252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9224077) q[2];
sx q[2];
rz(-1.2097404) q[2];
sx q[2];
rz(-2.6694471) q[2];
rz(1.2989429) q[3];
sx q[3];
rz(-1.8483714) q[3];
sx q[3];
rz(0.81056547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9206813) q[0];
sx q[0];
rz(-0.26316106) q[0];
sx q[0];
rz(-2.8780908) q[0];
rz(2.0384516) q[1];
sx q[1];
rz(-1.8270854) q[1];
sx q[1];
rz(-2.7679494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0804701) q[0];
sx q[0];
rz(-1.6310638) q[0];
sx q[0];
rz(-1.4979657) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2977826) q[2];
sx q[2];
rz(-1.4792974) q[2];
sx q[2];
rz(1.757427) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1800823) q[1];
sx q[1];
rz(-1.0133378) q[1];
sx q[1];
rz(2.7649759) q[1];
rz(2.3832537) q[3];
sx q[3];
rz(-0.74303526) q[3];
sx q[3];
rz(2.5287927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0282447) q[2];
sx q[2];
rz(-0.17013203) q[2];
sx q[2];
rz(2.587758) q[2];
rz(-1.3977741) q[3];
sx q[3];
rz(-0.60269409) q[3];
sx q[3];
rz(-2.8288614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58436191) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(1.9118017) q[0];
rz(2.9025485) q[1];
sx q[1];
rz(-1.6300647) q[1];
sx q[1];
rz(-2.8412433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2984788) q[0];
sx q[0];
rz(-1.4550147) q[0];
sx q[0];
rz(-0.16364574) q[0];
rz(0.70234583) q[2];
sx q[2];
rz(-1.1804198) q[2];
sx q[2];
rz(2.3927488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21619851) q[1];
sx q[1];
rz(-1.555519) q[1];
sx q[1];
rz(-3.141204) q[1];
x q[2];
rz(0.024552931) q[3];
sx q[3];
rz(-2.2835287) q[3];
sx q[3];
rz(-2.7803382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0099237) q[2];
sx q[2];
rz(-1.4574304) q[2];
sx q[2];
rz(2.8966676) q[2];
rz(-2.6217672) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(0.68827099) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35172611) q[0];
sx q[0];
rz(-0.34857294) q[0];
sx q[0];
rz(-1.9785731) q[0];
rz(-3.0746958) q[1];
sx q[1];
rz(-1.4935378) q[1];
sx q[1];
rz(1.012872) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6179919) q[0];
sx q[0];
rz(-0.80192425) q[0];
sx q[0];
rz(-2.2973934) q[0];
rz(-2.1383019) q[2];
sx q[2];
rz(-1.7219647) q[2];
sx q[2];
rz(-2.0625819) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.098274577) q[1];
sx q[1];
rz(-1.2407082) q[1];
sx q[1];
rz(-1.2374452) q[1];
rz(0.49533923) q[3];
sx q[3];
rz(-2.3916349) q[3];
sx q[3];
rz(-1.0487674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11671994) q[2];
sx q[2];
rz(-2.1634384) q[2];
sx q[2];
rz(0.32279521) q[2];
rz(-0.60574496) q[3];
sx q[3];
rz(-0.79234684) q[3];
sx q[3];
rz(2.7927223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.054319687) q[0];
sx q[0];
rz(-2.9948586) q[0];
sx q[0];
rz(3.1291381) q[0];
rz(0.74673486) q[1];
sx q[1];
rz(-0.92339271) q[1];
sx q[1];
rz(-2.8616203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0664205) q[0];
sx q[0];
rz(-0.1233347) q[0];
sx q[0];
rz(-1.146011) q[0];
rz(2.0048281) q[2];
sx q[2];
rz(-1.2723337) q[2];
sx q[2];
rz(3.0500183) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1713531) q[1];
sx q[1];
rz(-0.84988028) q[1];
sx q[1];
rz(1.9848787) q[1];
x q[2];
rz(0.2467201) q[3];
sx q[3];
rz(-1.9233875) q[3];
sx q[3];
rz(2.9143434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.81194699) q[2];
sx q[2];
rz(-2.3880366) q[2];
sx q[2];
rz(-2.9296056) q[2];
rz(0.82344615) q[3];
sx q[3];
rz(-1.4444838) q[3];
sx q[3];
rz(2.8823891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9987746) q[0];
sx q[0];
rz(-0.057567216) q[0];
sx q[0];
rz(0.69277358) q[0];
rz(-0.57299262) q[1];
sx q[1];
rz(-1.3447821) q[1];
sx q[1];
rz(-2.7105892) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7438223) q[0];
sx q[0];
rz(-1.5122688) q[0];
sx q[0];
rz(0.55214793) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0153322) q[2];
sx q[2];
rz(-1.017414) q[2];
sx q[2];
rz(-2.0972507) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.17934701) q[1];
sx q[1];
rz(-1.3192466) q[1];
sx q[1];
rz(0.72466447) q[1];
rz(1.0209084) q[3];
sx q[3];
rz(-1.3088994) q[3];
sx q[3];
rz(2.6481215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7196322) q[2];
sx q[2];
rz(-0.27406359) q[2];
sx q[2];
rz(-0.56023041) q[2];
rz(-2.6719921) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(-2.4277021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2847168) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(-0.84125413) q[1];
sx q[1];
rz(-2.0396736) q[1];
sx q[1];
rz(3.094818) q[1];
rz(-1.2011436) q[2];
sx q[2];
rz(-1.8664202) q[2];
sx q[2];
rz(2.0640434) q[2];
rz(0.91184323) q[3];
sx q[3];
rz(-1.5921436) q[3];
sx q[3];
rz(2.5719503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
