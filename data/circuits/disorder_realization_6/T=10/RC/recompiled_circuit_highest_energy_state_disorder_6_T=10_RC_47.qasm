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
rz(-2.1498635) q[0];
sx q[0];
rz(-2.5340134) q[0];
sx q[0];
rz(2.1460331) q[0];
rz(-2.6140656) q[1];
sx q[1];
rz(-1.0310643) q[1];
sx q[1];
rz(-0.48479015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3207389) q[0];
sx q[0];
rz(-1.0620097) q[0];
sx q[0];
rz(2.9910422) q[0];
rz(-pi) q[1];
rz(-0.44102168) q[2];
sx q[2];
rz(-0.82953757) q[2];
sx q[2];
rz(1.5462359) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.97580662) q[1];
sx q[1];
rz(-2.1790904) q[1];
sx q[1];
rz(3.0700141) q[1];
rz(2.9020314) q[3];
sx q[3];
rz(-2.8423841) q[3];
sx q[3];
rz(2.9394958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.506044) q[2];
sx q[2];
rz(-0.94749331) q[2];
sx q[2];
rz(0.49723899) q[2];
rz(-0.55073589) q[3];
sx q[3];
rz(-2.4516055) q[3];
sx q[3];
rz(1.1373038) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9841442) q[0];
sx q[0];
rz(-1.7406311) q[0];
sx q[0];
rz(-2.3781811) q[0];
rz(0.58452559) q[1];
sx q[1];
rz(-1.7598563) q[1];
sx q[1];
rz(2.1302628) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9824163) q[0];
sx q[0];
rz(-1.6859907) q[0];
sx q[0];
rz(-1.7027436) q[0];
rz(-pi) q[1];
rz(-1.6857851) q[2];
sx q[2];
rz(-1.3359112) q[2];
sx q[2];
rz(1.9180852) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8400795) q[1];
sx q[1];
rz(-0.40427819) q[1];
sx q[1];
rz(2.1188645) q[1];
x q[2];
rz(-0.39803466) q[3];
sx q[3];
rz(-2.0477437) q[3];
sx q[3];
rz(1.2881035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4804907) q[2];
sx q[2];
rz(-1.2031518) q[2];
sx q[2];
rz(2.1413546) q[2];
rz(0.58146042) q[3];
sx q[3];
rz(-2.4088819) q[3];
sx q[3];
rz(-1.940041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30705273) q[0];
sx q[0];
rz(-0.59891278) q[0];
sx q[0];
rz(0.55759984) q[0];
rz(-2.8343976) q[1];
sx q[1];
rz(-2.5418044) q[1];
sx q[1];
rz(0.44816005) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24009839) q[0];
sx q[0];
rz(-2.0848102) q[0];
sx q[0];
rz(-1.7230524) q[0];
rz(-pi) q[1];
rz(1.0501625) q[2];
sx q[2];
rz(-1.7562281) q[2];
sx q[2];
rz(0.38927573) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8413755) q[1];
sx q[1];
rz(-2.5102477) q[1];
sx q[1];
rz(1.6168827) q[1];
x q[2];
rz(0.43436111) q[3];
sx q[3];
rz(-1.54714) q[3];
sx q[3];
rz(-1.0756186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8223411) q[2];
sx q[2];
rz(-2.4288869) q[2];
sx q[2];
rz(2.5437497) q[2];
rz(-0.44761014) q[3];
sx q[3];
rz(-1.8748583) q[3];
sx q[3];
rz(-3.083057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7287801) q[0];
sx q[0];
rz(-3.1106115) q[0];
sx q[0];
rz(2.431751) q[0];
rz(0.9791044) q[1];
sx q[1];
rz(-2.282228) q[1];
sx q[1];
rz(-1.8202579) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13793547) q[0];
sx q[0];
rz(-0.70323479) q[0];
sx q[0];
rz(0.08589311) q[0];
x q[1];
rz(-0.15726451) q[2];
sx q[2];
rz(-2.3549821) q[2];
sx q[2];
rz(0.32341126) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.37934694) q[1];
sx q[1];
rz(-1.3031465) q[1];
sx q[1];
rz(-2.2257811) q[1];
rz(-pi) q[2];
rz(-2.0632486) q[3];
sx q[3];
rz(-0.60090099) q[3];
sx q[3];
rz(0.43481024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1506302) q[2];
sx q[2];
rz(-1.7143098) q[2];
sx q[2];
rz(2.4347351) q[2];
rz(-1.7957211) q[3];
sx q[3];
rz(-1.6790877) q[3];
sx q[3];
rz(0.70485392) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3020346) q[0];
sx q[0];
rz(-0.78080451) q[0];
sx q[0];
rz(1.6395521) q[0];
rz(1.543965) q[1];
sx q[1];
rz(-2.2815956) q[1];
sx q[1];
rz(1.4942716) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92163819) q[0];
sx q[0];
rz(-1.0003547) q[0];
sx q[0];
rz(0.0019689671) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59539184) q[2];
sx q[2];
rz(-1.5647581) q[2];
sx q[2];
rz(-1.572682) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8710821) q[1];
sx q[1];
rz(-1.8746261) q[1];
sx q[1];
rz(-0.31131502) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5389046) q[3];
sx q[3];
rz(-2.5261554) q[3];
sx q[3];
rz(-0.36356715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8813701) q[2];
sx q[2];
rz(-1.8689195) q[2];
sx q[2];
rz(0.55848813) q[2];
rz(2.6327366) q[3];
sx q[3];
rz(-1.5285834) q[3];
sx q[3];
rz(-1.8744899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5452071) q[0];
sx q[0];
rz(-1.243243) q[0];
sx q[0];
rz(0.16045706) q[0];
rz(2.4682553) q[1];
sx q[1];
rz(-1.2168177) q[1];
sx q[1];
rz(1.1268667) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25803963) q[0];
sx q[0];
rz(-2.6614688) q[0];
sx q[0];
rz(-0.86999805) q[0];
rz(-pi) q[1];
rz(2.8062988) q[2];
sx q[2];
rz(-0.79293434) q[2];
sx q[2];
rz(-1.060134) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1231251) q[1];
sx q[1];
rz(-1.3686124) q[1];
sx q[1];
rz(-2.1582161) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9301008) q[3];
sx q[3];
rz(-1.0387254) q[3];
sx q[3];
rz(-2.718975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1274073) q[2];
sx q[2];
rz(-1.5421474) q[2];
sx q[2];
rz(-2.3678153) q[2];
rz(-2.1524147) q[3];
sx q[3];
rz(-0.4072322) q[3];
sx q[3];
rz(-1.0729084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9147375) q[0];
sx q[0];
rz(-1.6586774) q[0];
sx q[0];
rz(0.7731272) q[0];
rz(1.5559366) q[1];
sx q[1];
rz(-1.8525886) q[1];
sx q[1];
rz(-1.9645436) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53922916) q[0];
sx q[0];
rz(-1.6022816) q[0];
sx q[0];
rz(-0.8336153) q[0];
rz(-0.99780166) q[2];
sx q[2];
rz(-1.3062451) q[2];
sx q[2];
rz(-0.71392347) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8985473) q[1];
sx q[1];
rz(-2.4345401) q[1];
sx q[1];
rz(-1.9999595) q[1];
rz(-pi) q[2];
rz(-0.1775976) q[3];
sx q[3];
rz(-1.4334205) q[3];
sx q[3];
rz(0.053973764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3561463) q[2];
sx q[2];
rz(-0.36474228) q[2];
sx q[2];
rz(2.6981603) q[2];
rz(-2.7555079) q[3];
sx q[3];
rz(-1.4177136) q[3];
sx q[3];
rz(-2.940787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76483738) q[0];
sx q[0];
rz(-1.5489464) q[0];
sx q[0];
rz(-3.1078597) q[0];
rz(-2.4434166) q[1];
sx q[1];
rz(-2.5444784) q[1];
sx q[1];
rz(0.66506344) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3515776) q[0];
sx q[0];
rz(-0.52886334) q[0];
sx q[0];
rz(1.0306503) q[0];
rz(-pi) q[1];
rz(-0.84455281) q[2];
sx q[2];
rz(-1.204927) q[2];
sx q[2];
rz(-1.6576115) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5387676) q[1];
sx q[1];
rz(-0.95315779) q[1];
sx q[1];
rz(-0.62299872) q[1];
rz(1.5499209) q[3];
sx q[3];
rz(-1.890278) q[3];
sx q[3];
rz(0.11939458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72208059) q[2];
sx q[2];
rz(-0.91700143) q[2];
sx q[2];
rz(-1.7783995) q[2];
rz(-2.6204387) q[3];
sx q[3];
rz(-1.8173953) q[3];
sx q[3];
rz(-0.49191973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1727961) q[0];
sx q[0];
rz(-1.7954614) q[0];
sx q[0];
rz(-2.3222493) q[0];
rz(0.68483812) q[1];
sx q[1];
rz(-1.2872773) q[1];
sx q[1];
rz(-3.0192764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4222833) q[0];
sx q[0];
rz(-0.78351142) q[0];
sx q[0];
rz(2.1859403) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1048466) q[2];
sx q[2];
rz(-2.4754562) q[2];
sx q[2];
rz(0.18830794) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.015339) q[1];
sx q[1];
rz(-2.1623134) q[1];
sx q[1];
rz(-1.9638318) q[1];
x q[2];
rz(-2.5531386) q[3];
sx q[3];
rz(-1.1752442) q[3];
sx q[3];
rz(1.5423397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.704432) q[2];
sx q[2];
rz(-1.3334583) q[2];
sx q[2];
rz(2.2740347) q[2];
rz(1.3982754) q[3];
sx q[3];
rz(-2.7531392) q[3];
sx q[3];
rz(-2.5625663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.687998) q[0];
sx q[0];
rz(-1.2247676) q[0];
sx q[0];
rz(1.0892185) q[0];
rz(-0.047529686) q[1];
sx q[1];
rz(-1.0462953) q[1];
sx q[1];
rz(-0.68738031) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9935038) q[0];
sx q[0];
rz(-1.5464791) q[0];
sx q[0];
rz(-2.0612102) q[0];
rz(-2.8609852) q[2];
sx q[2];
rz(-1.1550552) q[2];
sx q[2];
rz(3.0688697) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.43722758) q[1];
sx q[1];
rz(-0.95672551) q[1];
sx q[1];
rz(2.1350056) q[1];
rz(-pi) q[2];
x q[2];
rz(0.093133868) q[3];
sx q[3];
rz(-0.93619213) q[3];
sx q[3];
rz(2.840691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1462732) q[2];
sx q[2];
rz(-0.61890382) q[2];
sx q[2];
rz(1.6590365) q[2];
rz(2.4847374) q[3];
sx q[3];
rz(-1.1494136) q[3];
sx q[3];
rz(-0.38203865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8889846) q[0];
sx q[0];
rz(-1.1532619) q[0];
sx q[0];
rz(1.1070195) q[0];
rz(0.57869115) q[1];
sx q[1];
rz(-2.3042669) q[1];
sx q[1];
rz(-0.28490983) q[1];
rz(-2.9215521) q[2];
sx q[2];
rz(-2.7773428) q[2];
sx q[2];
rz(-0.11014948) q[2];
rz(-0.5169576) q[3];
sx q[3];
rz(-1.4993954) q[3];
sx q[3];
rz(-2.7837337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
