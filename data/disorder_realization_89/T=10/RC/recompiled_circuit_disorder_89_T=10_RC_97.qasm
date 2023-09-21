OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(-2.8470319) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.598782) q[0];
sx q[0];
rz(-3.071975) q[0];
sx q[0];
rz(1.9570062) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6099036) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(1.6959015) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1301073) q[1];
sx q[1];
rz(-1.3890146) q[1];
sx q[1];
rz(-0.45244042) q[1];
x q[2];
rz(0.59395091) q[3];
sx q[3];
rz(-1.4547326) q[3];
sx q[3];
rz(0.56967294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(-1.4953556) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(-1.4650311) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7709748) q[0];
sx q[0];
rz(-0.3807225) q[0];
sx q[0];
rz(-0.61852635) q[0];
x q[1];
rz(0.53441647) q[2];
sx q[2];
rz(-2.538531) q[2];
sx q[2];
rz(2.7433928) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6444781) q[1];
sx q[1];
rz(-1.2116355) q[1];
sx q[1];
rz(-0.70980806) q[1];
rz(-pi) q[2];
rz(-2.6133735) q[3];
sx q[3];
rz(-1.8209407) q[3];
sx q[3];
rz(1.6268829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(2.583288) q[2];
rz(1.0393418) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6168183) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(0.16361374) q[0];
rz(-0.71584654) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2046976) q[0];
sx q[0];
rz(-1.0442797) q[0];
sx q[0];
rz(1.599647) q[0];
rz(1.9203556) q[2];
sx q[2];
rz(-2.3974843) q[2];
sx q[2];
rz(-2.0958054) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.219017) q[1];
sx q[1];
rz(-2.1891928) q[1];
sx q[1];
rz(0.69505691) q[1];
rz(-pi) q[2];
rz(1.8009637) q[3];
sx q[3];
rz(-1.3535415) q[3];
sx q[3];
rz(2.1118856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(-1.019657) q[2];
rz(-0.9807469) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137988) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(0.83909488) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71826868) q[0];
sx q[0];
rz(-1.6800796) q[0];
sx q[0];
rz(-3.1302343) q[0];
x q[1];
rz(-2.7856366) q[2];
sx q[2];
rz(-1.9320556) q[2];
sx q[2];
rz(-1.661983) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8295633) q[1];
sx q[1];
rz(-2.2514572) q[1];
sx q[1];
rz(0.3373674) q[1];
rz(-pi) q[2];
rz(0.88704349) q[3];
sx q[3];
rz(-2.4409557) q[3];
sx q[3];
rz(0.62473434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.443976) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(0.237341) q[2];
rz(0.90732968) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.6326555) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(-1.6089815) q[0];
rz(-2.4081047) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(-1.3132494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35996138) q[0];
sx q[0];
rz(-0.93802035) q[0];
sx q[0];
rz(2.7347793) q[0];
x q[1];
rz(1.0079185) q[2];
sx q[2];
rz(-3.0745227) q[2];
sx q[2];
rz(-1.9474533) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5161799) q[1];
sx q[1];
rz(-1.1986294) q[1];
sx q[1];
rz(-2.7402997) q[1];
x q[2];
rz(-2.2458514) q[3];
sx q[3];
rz(-0.81479077) q[3];
sx q[3];
rz(-0.18273396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0017172) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(0.66429794) q[2];
rz(-0.70513606) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(-3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649081) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(0.054811906) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(2.6584113) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85855243) q[0];
sx q[0];
rz(-0.1917834) q[0];
sx q[0];
rz(-1.900308) q[0];
x q[1];
rz(-2.0334843) q[2];
sx q[2];
rz(-1.6001284) q[2];
sx q[2];
rz(-2.6372452) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.253787) q[1];
sx q[1];
rz(-2.0704381) q[1];
sx q[1];
rz(0.81412022) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9032848) q[3];
sx q[3];
rz(-1.4629435) q[3];
sx q[3];
rz(0.40963848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(-1.9412458) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(2.9912662) q[3];
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
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(0.95463395) q[0];
rz(-2.5700991) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(-0.25209299) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3207021) q[0];
sx q[0];
rz(-1.1517236) q[0];
sx q[0];
rz(1.7295895) q[0];
rz(2.3605395) q[2];
sx q[2];
rz(-2.1456246) q[2];
sx q[2];
rz(0.99036723) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8512307) q[1];
sx q[1];
rz(-0.92612672) q[1];
sx q[1];
rz(2.6840997) q[1];
x q[2];
rz(2.2456456) q[3];
sx q[3];
rz(-2.575867) q[3];
sx q[3];
rz(-1.5783527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(0.90448109) q[2];
rz(-1.0036184) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(1.7277539) q[0];
rz(0.66043234) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(-1.3716912) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28384128) q[0];
sx q[0];
rz(-2.5159266) q[0];
sx q[0];
rz(-1.7218504) q[0];
rz(-2.3805982) q[2];
sx q[2];
rz(-1.4174263) q[2];
sx q[2];
rz(-0.93190565) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.15118571) q[1];
sx q[1];
rz(-1.8670261) q[1];
sx q[1];
rz(2.5469261) q[1];
rz(-1.4189818) q[3];
sx q[3];
rz(-1.9278798) q[3];
sx q[3];
rz(-2.9861772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(0.39789847) q[2];
rz(-2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9234377) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(-1.6802457) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-2.871002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10588564) q[0];
sx q[0];
rz(-2.8286655) q[0];
sx q[0];
rz(-2.3653415) q[0];
rz(-pi) q[1];
rz(1.0762392) q[2];
sx q[2];
rz(-2.7521172) q[2];
sx q[2];
rz(1.4091834) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8447664) q[1];
sx q[1];
rz(-1.4335732) q[1];
sx q[1];
rz(1.574135) q[1];
rz(-3.0398265) q[3];
sx q[3];
rz(-1.093285) q[3];
sx q[3];
rz(-0.69672841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8018735) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(-1.6516997) q[2];
rz(-2.4387032) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-2.9872966) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-0.74238366) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61481793) q[0];
sx q[0];
rz(-0.97133884) q[0];
sx q[0];
rz(-2.4012698) q[0];
rz(-0.43138357) q[2];
sx q[2];
rz(-1.0215534) q[2];
sx q[2];
rz(-2.0768349) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87634516) q[1];
sx q[1];
rz(-2.7956388) q[1];
sx q[1];
rz(0.032851263) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5417682) q[3];
sx q[3];
rz(-1.4283984) q[3];
sx q[3];
rz(-0.27424973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2877038) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(1.3416946) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-1.2269536) q[2];
sx q[2];
rz(-0.77949066) q[2];
sx q[2];
rz(1.6463999) q[2];
rz(1.4349764) q[3];
sx q[3];
rz(-1.4321696) q[3];
sx q[3];
rz(-1.1925478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
