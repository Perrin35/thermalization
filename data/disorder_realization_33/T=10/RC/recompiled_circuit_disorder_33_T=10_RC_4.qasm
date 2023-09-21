OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(1.5703262) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9309064) q[0];
sx q[0];
rz(-1.9084198) q[0];
sx q[0];
rz(-2.7772285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3762796) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(-3.090976) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.624144) q[1];
sx q[1];
rz(-1.2341208) q[1];
sx q[1];
rz(1.8271853) q[1];
x q[2];
rz(2.9458463) q[3];
sx q[3];
rz(-1.2139075) q[3];
sx q[3];
rz(1.9407879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.87542614) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(1.1323294) q[2];
rz(1.4663565) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(-1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9448626) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(-0.18584132) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.8461684) q[1];
sx q[1];
rz(2.9247608) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8685535) q[0];
sx q[0];
rz(-2.2076063) q[0];
sx q[0];
rz(0.36006948) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88044135) q[2];
sx q[2];
rz(-2.464622) q[2];
sx q[2];
rz(2.0073236) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.20570457) q[1];
sx q[1];
rz(-1.7824031) q[1];
sx q[1];
rz(0.88797027) q[1];
rz(-pi) q[2];
rz(-2.5012245) q[3];
sx q[3];
rz(-2.159517) q[3];
sx q[3];
rz(2.0224188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8314787) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(1.2878093) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(1.3263946) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-0.93260971) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7639097) q[0];
sx q[0];
rz(-1.253486) q[0];
sx q[0];
rz(3.0850287) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7980174) q[2];
sx q[2];
rz(-1.7482687) q[2];
sx q[2];
rz(2.9858659) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4438666) q[1];
sx q[1];
rz(-1.6892471) q[1];
sx q[1];
rz(-0.5603793) q[1];
x q[2];
rz(2.1266537) q[3];
sx q[3];
rz(-1.1851289) q[3];
sx q[3];
rz(2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(2.0489342) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3595235) q[0];
sx q[0];
rz(-3.0451267) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(0.80530986) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(1.6436228) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1456137) q[0];
sx q[0];
rz(-1.0756452) q[0];
sx q[0];
rz(2.8061295) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0849886) q[2];
sx q[2];
rz(-1.7610234) q[2];
sx q[2];
rz(-0.20400001) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.30724635) q[1];
sx q[1];
rz(-1.3386968) q[1];
sx q[1];
rz(-0.26718617) q[1];
rz(-pi) q[2];
rz(-0.46521503) q[3];
sx q[3];
rz(-1.9808931) q[3];
sx q[3];
rz(-1.0614392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(0.70181075) q[2];
rz(-0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(-2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2410626) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(-0.81533122) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-1.048208) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74131504) q[0];
sx q[0];
rz(-1.321723) q[0];
sx q[0];
rz(-1.8427909) q[0];
rz(-1.5830718) q[2];
sx q[2];
rz(-0.91931146) q[2];
sx q[2];
rz(-2.2001681) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.198846) q[1];
sx q[1];
rz(-1.2578576) q[1];
sx q[1];
rz(-1.320977) q[1];
x q[2];
rz(-0.95543315) q[3];
sx q[3];
rz(-1.9121998) q[3];
sx q[3];
rz(-2.0379025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(-1.099951) q[2];
rz(2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0681756) q[0];
sx q[0];
rz(-0.59403479) q[0];
sx q[0];
rz(2.2391879) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-1.0598176) q[1];
sx q[1];
rz(-3.0117603) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3467305) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(2.5559588) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99545698) q[2];
sx q[2];
rz(-1.8838922) q[2];
sx q[2];
rz(0.42195937) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.063720238) q[1];
sx q[1];
rz(-1.2076326) q[1];
sx q[1];
rz(-0.6086463) q[1];
rz(-1.767166) q[3];
sx q[3];
rz(-1.0305627) q[3];
sx q[3];
rz(-1.2353209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-2.1924993) q[2];
sx q[2];
rz(-2.9373346) q[2];
rz(-1.2060818) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181353) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.4468505) q[0];
rz(-1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(-0.68626219) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2080363) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(-0.74525381) q[0];
rz(-pi) q[1];
rz(-2.3115736) q[2];
sx q[2];
rz(-1.1168715) q[2];
sx q[2];
rz(-1.1656851) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2346748) q[1];
sx q[1];
rz(-1.385681) q[1];
sx q[1];
rz(0.075637416) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5043082) q[3];
sx q[3];
rz(-2.6316959) q[3];
sx q[3];
rz(2.2989458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(3.1398204) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(-1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6034265) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.5015645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1743463) q[0];
sx q[0];
rz(-1.5902728) q[0];
sx q[0];
rz(1.0780225) q[0];
rz(-pi) q[1];
rz(1.6325475) q[2];
sx q[2];
rz(-1.5973063) q[2];
sx q[2];
rz(2.3960631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9403119) q[1];
sx q[1];
rz(-2.8021325) q[1];
sx q[1];
rz(-0.27075726) q[1];
rz(-pi) q[2];
rz(-2.4092259) q[3];
sx q[3];
rz(-2.3673956) q[3];
sx q[3];
rz(2.458651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(-1.2119279) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(0.31931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-2.7899172) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29031819) q[0];
sx q[0];
rz(-1.9418678) q[0];
sx q[0];
rz(-2.0593658) q[0];
rz(-pi) q[1];
rz(-2.4497689) q[2];
sx q[2];
rz(-1.8080538) q[2];
sx q[2];
rz(2.0123864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-2.3304686) q[1];
sx q[1];
rz(0.07304904) q[1];
rz(-pi) q[2];
rz(3.0319801) q[3];
sx q[3];
rz(-1.5188367) q[3];
sx q[3];
rz(0.51652858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7982771) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.6069501) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(1.0378029) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(2.1077572) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4170096) q[0];
sx q[0];
rz(-1.7010744) q[0];
sx q[0];
rz(2.2307322) q[0];
rz(2.1543703) q[2];
sx q[2];
rz(-2.3505031) q[2];
sx q[2];
rz(-0.9466048) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2227877) q[1];
sx q[1];
rz(-2.5606887) q[1];
sx q[1];
rz(1.7564299) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3528353) q[3];
sx q[3];
rz(-0.86943227) q[3];
sx q[3];
rz(-2.4234114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(0.6357843) q[2];
rz(-2.87129) q[3];
sx q[3];
rz(-2.342194) q[3];
sx q[3];
rz(-1.6132145) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6939659) q[0];
sx q[0];
rz(-1.3128558) q[0];
sx q[0];
rz(-2.0679612) q[0];
rz(-1.7059965) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-3.1096935) q[2];
sx q[2];
rz(-0.96822856) q[2];
sx q[2];
rz(-0.4005489) q[2];
rz(1.3676436) q[3];
sx q[3];
rz(-2.805134) q[3];
sx q[3];
rz(0.72611879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];