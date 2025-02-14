OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.91650668) q[0];
sx q[0];
rz(-2.9986311) q[0];
sx q[0];
rz(1.4885055) q[0];
rz(-2.0090964) q[1];
sx q[1];
rz(-0.43192616) q[1];
sx q[1];
rz(-2.5052524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36447517) q[0];
sx q[0];
rz(-2.416673) q[0];
sx q[0];
rz(0.97282797) q[0];
rz(-pi) q[1];
rz(-2.752334) q[2];
sx q[2];
rz(-1.2137652) q[2];
sx q[2];
rz(-0.87212002) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.072187034) q[1];
sx q[1];
rz(-1.9652307) q[1];
sx q[1];
rz(-0.13448127) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79145517) q[3];
sx q[3];
rz(-2.2558082) q[3];
sx q[3];
rz(0.58429694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0296313) q[2];
sx q[2];
rz(-1.1250857) q[2];
sx q[2];
rz(-0.65043989) q[2];
rz(1.316831) q[3];
sx q[3];
rz(-1.3464758) q[3];
sx q[3];
rz(-0.10183798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.0541075) q[0];
sx q[0];
rz(-0.58687812) q[0];
sx q[0];
rz(2.6869539) q[0];
rz(0.023160402) q[1];
sx q[1];
rz(-1.6975479) q[1];
sx q[1];
rz(-2.815411) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6070411) q[0];
sx q[0];
rz(-1.9605419) q[0];
sx q[0];
rz(2.5410065) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9942547) q[2];
sx q[2];
rz(-1.2988161) q[2];
sx q[2];
rz(-0.47951298) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4224976) q[1];
sx q[1];
rz(-1.0651673) q[1];
sx q[1];
rz(-2.0431594) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11740673) q[3];
sx q[3];
rz(-2.2281584) q[3];
sx q[3];
rz(-0.18169345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11639103) q[2];
sx q[2];
rz(-1.2025183) q[2];
sx q[2];
rz(2.1873059) q[2];
rz(1.8918234) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(-3.1402816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8702451) q[0];
sx q[0];
rz(-1.1529237) q[0];
sx q[0];
rz(1.9648319) q[0];
rz(-0.36508834) q[1];
sx q[1];
rz(-1.3026214) q[1];
sx q[1];
rz(-1.5286068) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919892) q[0];
sx q[0];
rz(-3.1168584) q[0];
sx q[0];
rz(-1.0036841) q[0];
x q[1];
rz(0.079054376) q[2];
sx q[2];
rz(-2.0131265) q[2];
sx q[2];
rz(2.4238264) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4903725) q[1];
sx q[1];
rz(-1.7224604) q[1];
sx q[1];
rz(-3.1107799) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5131475) q[3];
sx q[3];
rz(-1.639327) q[3];
sx q[3];
rz(1.203361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8927346) q[2];
sx q[2];
rz(-2.5881519) q[2];
sx q[2];
rz(1.152732) q[2];
rz(1.8966127) q[3];
sx q[3];
rz(-0.98074073) q[3];
sx q[3];
rz(-2.8583756) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1588441) q[0];
sx q[0];
rz(-11*pi/12) q[0];
sx q[0];
rz(-0.68956476) q[0];
rz(-3.116563) q[1];
sx q[1];
rz(-1.5182779) q[1];
sx q[1];
rz(1.7207346) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0293504) q[0];
sx q[0];
rz(-2.1922702) q[0];
sx q[0];
rz(-1.6981324) q[0];
rz(-pi) q[1];
rz(-2.9621012) q[2];
sx q[2];
rz(-1.9948261) q[2];
sx q[2];
rz(-1.9551203) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6384699) q[1];
sx q[1];
rz(-1.5984189) q[1];
sx q[1];
rz(-2.6383328) q[1];
rz(1.2989276) q[3];
sx q[3];
rz(-1.8966348) q[3];
sx q[3];
rz(-0.0001903521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.932852) q[2];
sx q[2];
rz(-2.9661861) q[2];
sx q[2];
rz(1.9746732) q[2];
rz(2.6973727) q[3];
sx q[3];
rz(-1.2362213) q[3];
sx q[3];
rz(0.3096295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.8532448) q[0];
sx q[0];
rz(-1.3695559) q[0];
sx q[0];
rz(0.41611588) q[0];
rz(-2.6314645) q[1];
sx q[1];
rz(-0.77113873) q[1];
sx q[1];
rz(-2.3663734) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907523) q[0];
sx q[0];
rz(-1.3715944) q[0];
sx q[0];
rz(0.61934031) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7367896) q[2];
sx q[2];
rz(-2.7049613) q[2];
sx q[2];
rz(-0.2704119) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9265799) q[1];
sx q[1];
rz(-1.6832507) q[1];
sx q[1];
rz(1.4447029) q[1];
rz(0.59501641) q[3];
sx q[3];
rz(-1.887768) q[3];
sx q[3];
rz(1.1770615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.2148718) q[2];
sx q[2];
rz(-0.99433172) q[2];
sx q[2];
rz(-2.1057687) q[2];
rz(2.9521247) q[3];
sx q[3];
rz(-0.47179705) q[3];
sx q[3];
rz(2.6389879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(1.254461) q[0];
sx q[0];
rz(-3.0939565) q[0];
sx q[0];
rz(0.3558085) q[0];
rz(1.4954781) q[1];
sx q[1];
rz(-2.1864086) q[1];
sx q[1];
rz(1.1411508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0063067181) q[0];
sx q[0];
rz(-1.7809488) q[0];
sx q[0];
rz(-2.2457613) q[0];
x q[1];
rz(3.1325794) q[2];
sx q[2];
rz(-2.1823332) q[2];
sx q[2];
rz(-1.0338986) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.49887744) q[1];
sx q[1];
rz(-1.6756454) q[1];
sx q[1];
rz(0.14123209) q[1];
rz(-1.5445956) q[3];
sx q[3];
rz(-2.7075534) q[3];
sx q[3];
rz(1.7569913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.63153875) q[2];
sx q[2];
rz(-2.4271991) q[2];
sx q[2];
rz(1.4368524) q[2];
rz(1.0531744) q[3];
sx q[3];
rz(-1.6097693) q[3];
sx q[3];
rz(-0.50301445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18008867) q[0];
sx q[0];
rz(-2.8453974) q[0];
sx q[0];
rz(0.12251138) q[0];
rz(2.4854614) q[1];
sx q[1];
rz(-1.8729112) q[1];
sx q[1];
rz(-1.3414541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2264474) q[0];
sx q[0];
rz(-1.6065803) q[0];
sx q[0];
rz(3.0962178) q[0];
rz(-1.4627005) q[2];
sx q[2];
rz(-1.8524395) q[2];
sx q[2];
rz(0.070527129) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9751898) q[1];
sx q[1];
rz(-2.4619048) q[1];
sx q[1];
rz(2.0026592) q[1];
rz(-pi) q[2];
rz(2.6977067) q[3];
sx q[3];
rz(-1.4680913) q[3];
sx q[3];
rz(0.041180276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.573367) q[2];
sx q[2];
rz(-1.2218916) q[2];
sx q[2];
rz(-1.4637671) q[2];
rz(2.407414) q[3];
sx q[3];
rz(-2.8575183) q[3];
sx q[3];
rz(-1.3045788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463335) q[0];
sx q[0];
rz(-1.7743552) q[0];
sx q[0];
rz(-0.45860589) q[0];
rz(-1.0147702) q[1];
sx q[1];
rz(-1.9472803) q[1];
sx q[1];
rz(1.0626622) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.053772702) q[0];
sx q[0];
rz(-1.7369978) q[0];
sx q[0];
rz(-1.6292162) q[0];
x q[1];
rz(-2.1193211) q[2];
sx q[2];
rz(-2.3000225) q[2];
sx q[2];
rz(2.4265576) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52673756) q[1];
sx q[1];
rz(-1.4397845) q[1];
sx q[1];
rz(3.0796771) q[1];
rz(0.32286947) q[3];
sx q[3];
rz(-2.0677276) q[3];
sx q[3];
rz(2.5931155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8899272) q[2];
sx q[2];
rz(-1.3883611) q[2];
sx q[2];
rz(1.3512705) q[2];
rz(-1.9225559) q[3];
sx q[3];
rz(-0.42870298) q[3];
sx q[3];
rz(-2.3082699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625921) q[0];
sx q[0];
rz(-1.3840249) q[0];
sx q[0];
rz(0.90743995) q[0];
rz(-0.22706789) q[1];
sx q[1];
rz(-1.0457958) q[1];
sx q[1];
rz(-2.2423832) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9010021) q[0];
sx q[0];
rz(-2.0792897) q[0];
sx q[0];
rz(2.2700538) q[0];
rz(0.078446968) q[2];
sx q[2];
rz(-1.4885934) q[2];
sx q[2];
rz(-1.3644827) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8912298) q[1];
sx q[1];
rz(-2.0768407) q[1];
sx q[1];
rz(0.21234588) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.090631) q[3];
sx q[3];
rz(-0.82333857) q[3];
sx q[3];
rz(2.2281102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1850618) q[2];
sx q[2];
rz(-1.1659634) q[2];
sx q[2];
rz(-2.5409307) q[2];
rz(1.0968084) q[3];
sx q[3];
rz(-2.5513702) q[3];
sx q[3];
rz(0.87305951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3908865) q[0];
sx q[0];
rz(-1.0989256) q[0];
sx q[0];
rz(-0.98439687) q[0];
rz(-0.4920494) q[1];
sx q[1];
rz(-1.4130519) q[1];
sx q[1];
rz(1.1955998) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0988783) q[0];
sx q[0];
rz(-0.98473583) q[0];
sx q[0];
rz(-1.1436966) q[0];
x q[1];
rz(-2.6457328) q[2];
sx q[2];
rz(-0.95520077) q[2];
sx q[2];
rz(0.95814182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1164536) q[1];
sx q[1];
rz(-2.7500116) q[1];
sx q[1];
rz(1.3571597) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77287425) q[3];
sx q[3];
rz(-1.9108655) q[3];
sx q[3];
rz(3.056385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0857346) q[2];
sx q[2];
rz(-1.6207638) q[2];
sx q[2];
rz(-0.98304191) q[2];
rz(-0.25162697) q[3];
sx q[3];
rz(-0.42767891) q[3];
sx q[3];
rz(1.0035286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6428103) q[0];
sx q[0];
rz(-2.1693873) q[0];
sx q[0];
rz(0.65709773) q[0];
rz(1.2596754) q[1];
sx q[1];
rz(-1.4114264) q[1];
sx q[1];
rz(-2.2687601) q[1];
rz(-0.85763422) q[2];
sx q[2];
rz(-1.7339755) q[2];
sx q[2];
rz(-0.82790464) q[2];
rz(0.76446492) q[3];
sx q[3];
rz(-0.44776147) q[3];
sx q[3];
rz(0.54816435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
