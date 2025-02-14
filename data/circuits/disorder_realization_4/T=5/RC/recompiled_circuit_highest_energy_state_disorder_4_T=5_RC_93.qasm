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
rz(0.29851222) q[0];
sx q[0];
rz(-1.8520344) q[0];
sx q[0];
rz(0.02331743) q[0];
rz(-2.159637) q[1];
sx q[1];
rz(-2.4040931) q[1];
sx q[1];
rz(1.6657383) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97853547) q[0];
sx q[0];
rz(-0.082162372) q[0];
sx q[0];
rz(-2.1294247) q[0];
rz(-pi) q[1];
rz(-0.040272399) q[2];
sx q[2];
rz(-1.9007517) q[2];
sx q[2];
rz(-1.2023991) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9495342) q[1];
sx q[1];
rz(-0.84334521) q[1];
sx q[1];
rz(0.53967584) q[1];
x q[2];
rz(-0.9528927) q[3];
sx q[3];
rz(-0.77761071) q[3];
sx q[3];
rz(-0.67985204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4618571) q[2];
sx q[2];
rz(-1.4718141) q[2];
sx q[2];
rz(-0.30065817) q[2];
rz(2.4833208) q[3];
sx q[3];
rz(-0.39977795) q[3];
sx q[3];
rz(2.5532653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41475007) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(-0.17587371) q[0];
rz(-2.580592) q[1];
sx q[1];
rz(-2.439552) q[1];
sx q[1];
rz(2.9452513) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4787653) q[0];
sx q[0];
rz(-1.9050373) q[0];
sx q[0];
rz(-1.1007376) q[0];
rz(-0.11949338) q[2];
sx q[2];
rz(-1.4548848) q[2];
sx q[2];
rz(1.8279861) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.92561103) q[1];
sx q[1];
rz(-2.6261733) q[1];
sx q[1];
rz(-1.5748596) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7815468) q[3];
sx q[3];
rz(-2.776675) q[3];
sx q[3];
rz(1.5555738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.089513) q[2];
sx q[2];
rz(-0.61605805) q[2];
sx q[2];
rz(-1.4698131) q[2];
rz(2.6164264) q[3];
sx q[3];
rz(-1.4273806) q[3];
sx q[3];
rz(-2.4770881) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0199652) q[0];
sx q[0];
rz(-0.88931924) q[0];
sx q[0];
rz(2.5110733) q[0];
rz(-1.1446674) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(-3.0847881) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5057748) q[0];
sx q[0];
rz(-2.5774101) q[0];
sx q[0];
rz(2.3057372) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1326829) q[2];
sx q[2];
rz(-0.76912921) q[2];
sx q[2];
rz(3.1388856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.60384974) q[1];
sx q[1];
rz(-1.7122388) q[1];
sx q[1];
rz(0.46700041) q[1];
x q[2];
rz(-1.9980691) q[3];
sx q[3];
rz(-1.7649073) q[3];
sx q[3];
rz(2.8097092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8593665) q[2];
sx q[2];
rz(-1.4713918) q[2];
sx q[2];
rz(2.7024506) q[2];
rz(-1.3965083) q[3];
sx q[3];
rz(-1.2625932) q[3];
sx q[3];
rz(2.5631574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0466995) q[0];
sx q[0];
rz(-2.4201604) q[0];
sx q[0];
rz(-1.0149581) q[0];
rz(3.0789442) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(0.28422022) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5361523) q[0];
sx q[0];
rz(-1.9198737) q[0];
sx q[0];
rz(-0.45850584) q[0];
rz(-pi) q[1];
rz(1.7106871) q[2];
sx q[2];
rz(-2.0645185) q[2];
sx q[2];
rz(1.3700652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9584459) q[1];
sx q[1];
rz(-1.5827521) q[1];
sx q[1];
rz(-2.7832116) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3174414) q[3];
sx q[3];
rz(-1.5572963) q[3];
sx q[3];
rz(-0.48397348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3011498) q[2];
sx q[2];
rz(-1.0089077) q[2];
sx q[2];
rz(-0.76048771) q[2];
rz(2.8216951) q[3];
sx q[3];
rz(-2.0127998) q[3];
sx q[3];
rz(2.1130051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29192057) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(-0.40529761) q[0];
rz(-0.48794508) q[1];
sx q[1];
rz(-1.1064203) q[1];
sx q[1];
rz(2.778756) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80984801) q[0];
sx q[0];
rz(-0.54143006) q[0];
sx q[0];
rz(-2.5141513) q[0];
x q[1];
rz(-1.9486729) q[2];
sx q[2];
rz(-2.123017) q[2];
sx q[2];
rz(3.0966126) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88089534) q[1];
sx q[1];
rz(-2.915307) q[1];
sx q[1];
rz(-2.923171) q[1];
x q[2];
rz(1.8764349) q[3];
sx q[3];
rz(-0.89857093) q[3];
sx q[3];
rz(1.6068939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1399347) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(-2.4190767) q[2];
rz(0.51269382) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(1.9374013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2672511) q[0];
sx q[0];
rz(-2.5047472) q[0];
sx q[0];
rz(-0.27477086) q[0];
rz(0.35708669) q[1];
sx q[1];
rz(-1.6386702) q[1];
sx q[1];
rz(-1.9836609) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0561075) q[0];
sx q[0];
rz(-0.72626136) q[0];
sx q[0];
rz(1.1994491) q[0];
x q[1];
rz(-2.2311796) q[2];
sx q[2];
rz(-1.0779625) q[2];
sx q[2];
rz(-3.1255258) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5722583) q[1];
sx q[1];
rz(-0.63779325) q[1];
sx q[1];
rz(0.34697726) q[1];
x q[2];
rz(0.66863184) q[3];
sx q[3];
rz(-0.9210081) q[3];
sx q[3];
rz(-1.1710407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.57924119) q[2];
sx q[2];
rz(-1.7513821) q[2];
sx q[2];
rz(-0.58970279) q[2];
rz(1.3127182) q[3];
sx q[3];
rz(-1.6629985) q[3];
sx q[3];
rz(2.5780799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8967459) q[0];
sx q[0];
rz(-2.1488996) q[0];
sx q[0];
rz(0.27031159) q[0];
rz(2.7145794) q[1];
sx q[1];
rz(-1.3753563) q[1];
sx q[1];
rz(-0.40831533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6231096) q[0];
sx q[0];
rz(-0.38869959) q[0];
sx q[0];
rz(-2.8220909) q[0];
rz(-0.87323453) q[2];
sx q[2];
rz(-1.6484652) q[2];
sx q[2];
rz(-0.49688646) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6417676) q[1];
sx q[1];
rz(-2.1534792) q[1];
sx q[1];
rz(-0.50104435) q[1];
rz(-pi) q[2];
rz(-0.116577) q[3];
sx q[3];
rz(-0.72598476) q[3];
sx q[3];
rz(0.67949142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.938505) q[2];
sx q[2];
rz(-1.0838584) q[2];
sx q[2];
rz(2.2173524) q[2];
rz(0.83485323) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(0.99669641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1860109) q[0];
sx q[0];
rz(-0.38014933) q[0];
sx q[0];
rz(-2.7741449) q[0];
rz(0.5683178) q[1];
sx q[1];
rz(-1.3325997) q[1];
sx q[1];
rz(-0.41499358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58632219) q[0];
sx q[0];
rz(-1.8316226) q[0];
sx q[0];
rz(-2.0404979) q[0];
x q[1];
rz(-1.0239564) q[2];
sx q[2];
rz(-1.1073565) q[2];
sx q[2];
rz(-1.6187504) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6549544) q[1];
sx q[1];
rz(-2.1732501) q[1];
sx q[1];
rz(2.2008373) q[1];
rz(-pi) q[2];
rz(1.9165809) q[3];
sx q[3];
rz(-1.5756851) q[3];
sx q[3];
rz(-1.1556311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35855287) q[2];
sx q[2];
rz(-2.2542605) q[2];
sx q[2];
rz(0.83317327) q[2];
rz(2.7821275) q[3];
sx q[3];
rz(-2.0514252) q[3];
sx q[3];
rz(1.5045213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1737162) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(0.014884431) q[0];
rz(1.124565) q[1];
sx q[1];
rz(-0.24528565) q[1];
sx q[1];
rz(-3.0344149) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10725645) q[0];
sx q[0];
rz(-1.7198623) q[0];
sx q[0];
rz(2.3019522) q[0];
x q[1];
rz(-2.066266) q[2];
sx q[2];
rz(-2.3297254) q[2];
sx q[2];
rz(-0.56174413) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31871492) q[1];
sx q[1];
rz(-1.7194304) q[1];
sx q[1];
rz(-2.6591402) q[1];
rz(-pi) q[2];
rz(1.393287) q[3];
sx q[3];
rz(-1.313398) q[3];
sx q[3];
rz(2.7696115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5643481) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(-2.2819819) q[2];
rz(0.25935069) q[3];
sx q[3];
rz(-1.3602463) q[3];
sx q[3];
rz(0.41659659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49941007) q[0];
sx q[0];
rz(-3.0152617) q[0];
sx q[0];
rz(-1.9835749) q[0];
rz(0.49531373) q[1];
sx q[1];
rz(-1.1664349) q[1];
sx q[1];
rz(0.29702979) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3756145) q[0];
sx q[0];
rz(-2.8225464) q[0];
sx q[0];
rz(-1.5837529) q[0];
rz(1.6112686) q[2];
sx q[2];
rz(-1.2707316) q[2];
sx q[2];
rz(-1.435868) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2386855) q[1];
sx q[1];
rz(-2.8082962) q[1];
sx q[1];
rz(-2.6413647) q[1];
rz(3.1252485) q[3];
sx q[3];
rz(-0.90829231) q[3];
sx q[3];
rz(-1.3682883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.26934066) q[2];
sx q[2];
rz(-0.38186914) q[2];
sx q[2];
rz(-0.86088172) q[2];
rz(2.6779029) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(2.004682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9853482) q[0];
sx q[0];
rz(-1.5679659) q[0];
sx q[0];
rz(1.1563942) q[0];
rz(-1.8078177) q[1];
sx q[1];
rz(-1.6009686) q[1];
sx q[1];
rz(0.70645465) q[1];
rz(0.68177022) q[2];
sx q[2];
rz(-2.3538156) q[2];
sx q[2];
rz(0.71629477) q[2];
rz(1.1838589) q[3];
sx q[3];
rz(-1.1087742) q[3];
sx q[3];
rz(-2.5840989) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
