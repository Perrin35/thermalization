OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7528485) q[0];
sx q[0];
rz(-0.53628439) q[0];
sx q[0];
rz(2.1938238) q[0];
rz(-1.3287969) q[1];
sx q[1];
rz(-1.8741908) q[1];
sx q[1];
rz(1.0277494) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2078903) q[0];
sx q[0];
rz(-2.7650802) q[0];
sx q[0];
rz(-0.062113751) q[0];
x q[1];
rz(-2.920354) q[2];
sx q[2];
rz(-2.049963) q[2];
sx q[2];
rz(-2.0146807) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6537522) q[1];
sx q[1];
rz(-2.2334705) q[1];
sx q[1];
rz(-2.0648048) q[1];
x q[2];
rz(0.62698934) q[3];
sx q[3];
rz(-1.1999745) q[3];
sx q[3];
rz(-0.8644608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(-1.1132647) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(0.072409078) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(-2.8438399) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(-2.0334977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0034804) q[0];
sx q[0];
rz(-1.7797884) q[0];
sx q[0];
rz(-0.94603993) q[0];
rz(2.1531395) q[2];
sx q[2];
rz(-2.4465912) q[2];
sx q[2];
rz(0.37441355) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.019921692) q[1];
sx q[1];
rz(-2.6056075) q[1];
sx q[1];
rz(-2.3977445) q[1];
rz(-1.3012582) q[3];
sx q[3];
rz(-2.2552935) q[3];
sx q[3];
rz(1.0052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(0.97529808) q[2];
rz(2.2235928) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(-2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(-2.1767445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7643825) q[0];
sx q[0];
rz(-1.2493734) q[0];
sx q[0];
rz(-1.8660603) q[0];
x q[1];
rz(1.6134878) q[2];
sx q[2];
rz(-1.2741158) q[2];
sx q[2];
rz(1.3670849) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.94273401) q[1];
sx q[1];
rz(-0.13730362) q[1];
sx q[1];
rz(1.0820461) q[1];
rz(-pi) q[2];
rz(-0.19823719) q[3];
sx q[3];
rz(-1.9803515) q[3];
sx q[3];
rz(-0.73892456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(2.9099693) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(-2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(2.8994697) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0771368) q[0];
sx q[0];
rz(-2.4479439) q[0];
sx q[0];
rz(2.8855188) q[0];
x q[1];
rz(-0.42963117) q[2];
sx q[2];
rz(-1.5589082) q[2];
sx q[2];
rz(2.5846543) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.081269216) q[1];
sx q[1];
rz(-2.6594866) q[1];
sx q[1];
rz(2.768885) q[1];
rz(-pi) q[2];
rz(2.8321213) q[3];
sx q[3];
rz(-0.43635338) q[3];
sx q[3];
rz(2.5040124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1233998) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-3.0467395) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(2.7807996) q[0];
rz(1.3882673) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(1.1345908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.40588) q[0];
sx q[0];
rz(-0.94731936) q[0];
sx q[0];
rz(0.9491802) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0005433) q[2];
sx q[2];
rz(-1.8998713) q[2];
sx q[2];
rz(2.3171901) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7736588) q[1];
sx q[1];
rz(-1.415167) q[1];
sx q[1];
rz(-2.6917798) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9633425) q[3];
sx q[3];
rz(-0.44294391) q[3];
sx q[3];
rz(-0.49803842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0052884) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(2.5689382) q[2];
rz(-0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3271493) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(1.2231187) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(1.1522326) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5778225) q[0];
sx q[0];
rz(-1.676079) q[0];
sx q[0];
rz(3.1288414) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67201519) q[2];
sx q[2];
rz(-2.5294371) q[2];
sx q[2];
rz(-1.9415346) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7631543) q[1];
sx q[1];
rz(-1.7708781) q[1];
sx q[1];
rz(1.6775908) q[1];
rz(-pi) q[2];
rz(-0.96529393) q[3];
sx q[3];
rz(-0.98635841) q[3];
sx q[3];
rz(2.8402929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46889177) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(3.0793072) q[0];
rz(0.1858055) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(2.7468162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4515526) q[0];
sx q[0];
rz(-1.3344904) q[0];
sx q[0];
rz(2.7355746) q[0];
x q[1];
rz(-2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(-1.2303908) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2784087) q[1];
sx q[1];
rz(-1.8179968) q[1];
sx q[1];
rz(2.7017038) q[1];
rz(-1.5338321) q[3];
sx q[3];
rz(-1.7241524) q[3];
sx q[3];
rz(-1.6085898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(-0.27302343) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(2.8295512) q[3];
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
rz(-2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(-1.2851108) q[0];
rz(1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(-1.3407019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1788951) q[0];
sx q[0];
rz(-0.51998752) q[0];
sx q[0];
rz(0.87332256) q[0];
rz(-pi) q[1];
x q[1];
rz(0.060872002) q[2];
sx q[2];
rz(-2.4448622) q[2];
sx q[2];
rz(-0.75887647) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.50780523) q[1];
sx q[1];
rz(-1.2101189) q[1];
sx q[1];
rz(-2.7651869) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1741381) q[3];
sx q[3];
rz(-0.52166044) q[3];
sx q[3];
rz(-2.012901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(1.0236053) q[2];
rz(0.18493955) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.0969365) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(1.6212844) q[0];
rz(0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(2.3044589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1820113) q[0];
sx q[0];
rz(-0.89996979) q[0];
sx q[0];
rz(1.8295893) q[0];
rz(-pi) q[1];
rz(0.43948549) q[2];
sx q[2];
rz(-0.62289933) q[2];
sx q[2];
rz(-0.57657951) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76444641) q[1];
sx q[1];
rz(-0.18491491) q[1];
sx q[1];
rz(-0.96692337) q[1];
rz(2.0765452) q[3];
sx q[3];
rz(-2.5878083) q[3];
sx q[3];
rz(-1.5050642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5403486) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(-0.17962757) q[2];
rz(2.1458697) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(-1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(-2.0843704) q[0];
rz(-3.0341042) q[1];
sx q[1];
rz(-1.8881533) q[1];
sx q[1];
rz(2.1616139) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5532593) q[0];
sx q[0];
rz(-0.13561121) q[0];
sx q[0];
rz(1.9562264) q[0];
rz(-pi) q[1];
rz(2.6238742) q[2];
sx q[2];
rz(-1.1283518) q[2];
sx q[2];
rz(0.25416086) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6085538) q[1];
sx q[1];
rz(-1.8095008) q[1];
sx q[1];
rz(2.0134316) q[1];
rz(-pi) q[2];
rz(1.8633217) q[3];
sx q[3];
rz(-2.8046126) q[3];
sx q[3];
rz(-1.8388336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3048627) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(1.0277964) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733611) q[0];
sx q[0];
rz(-1.0366476) q[0];
sx q[0];
rz(-1.114053) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-0.17765799) q[2];
sx q[2];
rz(-1.0189609) q[2];
sx q[2];
rz(-0.56448274) q[2];
rz(2.6408623) q[3];
sx q[3];
rz(-1.9095608) q[3];
sx q[3];
rz(0.6469971) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
