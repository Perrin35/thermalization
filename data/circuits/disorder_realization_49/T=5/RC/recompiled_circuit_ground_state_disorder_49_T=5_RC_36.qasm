OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9743118) q[0];
sx q[0];
rz(-0.33742961) q[0];
sx q[0];
rz(0.20198241) q[0];
rz(2.1972411) q[1];
sx q[1];
rz(-0.75466067) q[1];
sx q[1];
rz(1.4442297) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.429753) q[0];
sx q[0];
rz(-2.0637207) q[0];
sx q[0];
rz(-1.3322387) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0488705) q[2];
sx q[2];
rz(-0.23140027) q[2];
sx q[2];
rz(-2.8211942) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66975245) q[1];
sx q[1];
rz(-0.83263157) q[1];
sx q[1];
rz(0.91186146) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0962288) q[3];
sx q[3];
rz(-1.9352018) q[3];
sx q[3];
rz(1.7045329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.28295383) q[2];
sx q[2];
rz(-0.084246548) q[2];
sx q[2];
rz(1.5343182) q[2];
rz(0.18680799) q[3];
sx q[3];
rz(-2.6818633) q[3];
sx q[3];
rz(1.4927347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9635791) q[0];
sx q[0];
rz(-2.3602965) q[0];
sx q[0];
rz(-0.70469967) q[0];
rz(2.1251382) q[1];
sx q[1];
rz(-1.1926788) q[1];
sx q[1];
rz(1.1526398) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4999324) q[0];
sx q[0];
rz(-1.0295273) q[0];
sx q[0];
rz(-0.199078) q[0];
rz(-2.1179885) q[2];
sx q[2];
rz(-2.0788801) q[2];
sx q[2];
rz(1.91712) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.94344596) q[1];
sx q[1];
rz(-1.0662765) q[1];
sx q[1];
rz(2.2322502) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14039881) q[3];
sx q[3];
rz(-0.75274668) q[3];
sx q[3];
rz(1.5206159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.54027259) q[2];
sx q[2];
rz(-2.2728964) q[2];
sx q[2];
rz(1.3199838) q[2];
rz(0.071062239) q[3];
sx q[3];
rz(-0.080852121) q[3];
sx q[3];
rz(-1.094187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1048626) q[0];
sx q[0];
rz(-2.7920089) q[0];
sx q[0];
rz(2.2678251) q[0];
rz(-1.8050487) q[1];
sx q[1];
rz(-2.4599894) q[1];
sx q[1];
rz(2.2866586) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.076043) q[0];
sx q[0];
rz(-2.4618755) q[0];
sx q[0];
rz(1.9125841) q[0];
x q[1];
rz(-2.1941691) q[2];
sx q[2];
rz(-0.50259841) q[2];
sx q[2];
rz(2.0386774) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1347772) q[1];
sx q[1];
rz(-1.2689277) q[1];
sx q[1];
rz(1.6848438) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.072238) q[3];
sx q[3];
rz(-2.4756458) q[3];
sx q[3];
rz(-1.7518821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9285589) q[2];
sx q[2];
rz(-2.010767) q[2];
sx q[2];
rz(2.0862759) q[2];
rz(-0.94163752) q[3];
sx q[3];
rz(-1.0712653) q[3];
sx q[3];
rz(2.1987703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4554491) q[0];
sx q[0];
rz(-0.84857714) q[0];
sx q[0];
rz(0.69547478) q[0];
rz(1.4675568) q[1];
sx q[1];
rz(-1.8753884) q[1];
sx q[1];
rz(3.0470336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7692663) q[0];
sx q[0];
rz(-2.8594198) q[0];
sx q[0];
rz(1.7442987) q[0];
x q[1];
rz(0.60463011) q[2];
sx q[2];
rz(-2.1994091) q[2];
sx q[2];
rz(0.18456799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.27457224) q[1];
sx q[1];
rz(-0.34718514) q[1];
sx q[1];
rz(-0.95496655) q[1];
rz(-pi) q[2];
rz(-2.2587682) q[3];
sx q[3];
rz(-1.7738288) q[3];
sx q[3];
rz(-2.1830851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5203984) q[2];
sx q[2];
rz(-0.79402557) q[2];
sx q[2];
rz(-2.3801079) q[2];
rz(2.9705808) q[3];
sx q[3];
rz(-1.8595502) q[3];
sx q[3];
rz(-3.0448992) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3068202) q[0];
sx q[0];
rz(-0.67245317) q[0];
sx q[0];
rz(2.6493454) q[0];
rz(-1.3193839) q[1];
sx q[1];
rz(-1.4232114) q[1];
sx q[1];
rz(0.48670235) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4159631) q[0];
sx q[0];
rz(-1.6902903) q[0];
sx q[0];
rz(2.7404345) q[0];
rz(-2.7560744) q[2];
sx q[2];
rz(-0.80759128) q[2];
sx q[2];
rz(-2.453192) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.25025955) q[1];
sx q[1];
rz(-0.57602611) q[1];
sx q[1];
rz(-0.86153077) q[1];
rz(-1.1718122) q[3];
sx q[3];
rz(-0.55087844) q[3];
sx q[3];
rz(-0.18607947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95185602) q[2];
sx q[2];
rz(-2.2978013) q[2];
sx q[2];
rz(1.0978511) q[2];
rz(-0.83419937) q[3];
sx q[3];
rz(-0.67535496) q[3];
sx q[3];
rz(1.437291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.6247691) q[0];
sx q[0];
rz(-0.56684816) q[0];
sx q[0];
rz(0.045850642) q[0];
rz(-2.397873) q[1];
sx q[1];
rz(-0.36879483) q[1];
sx q[1];
rz(-3.0812982) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9550388) q[0];
sx q[0];
rz(-2.5738724) q[0];
sx q[0];
rz(-2.4283152) q[0];
rz(-pi) q[1];
rz(0.8745114) q[2];
sx q[2];
rz(-1.881045) q[2];
sx q[2];
rz(-0.53521148) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0963147) q[1];
sx q[1];
rz(-1.244925) q[1];
sx q[1];
rz(-2.9067944) q[1];
x q[2];
rz(0.89222091) q[3];
sx q[3];
rz(-1.4961317) q[3];
sx q[3];
rz(-1.1589662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9149949) q[2];
sx q[2];
rz(-1.8151585) q[2];
sx q[2];
rz(-2.162852) q[2];
rz(1.2616875) q[3];
sx q[3];
rz(-1.0190957) q[3];
sx q[3];
rz(-0.88920465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9642445) q[0];
sx q[0];
rz(-1.6151936) q[0];
sx q[0];
rz(0.88441315) q[0];
rz(-2.2824967) q[1];
sx q[1];
rz(-1.1235378) q[1];
sx q[1];
rz(-2.9335847) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7743083) q[0];
sx q[0];
rz(-1.2173664) q[0];
sx q[0];
rz(-2.6103781) q[0];
rz(2.5701373) q[2];
sx q[2];
rz(-0.68585912) q[2];
sx q[2];
rz(-1.4446007) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.35321843) q[1];
sx q[1];
rz(-0.93868449) q[1];
sx q[1];
rz(-0.33003426) q[1];
x q[2];
rz(0.0020387928) q[3];
sx q[3];
rz(-1.0674132) q[3];
sx q[3];
rz(0.17755213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7615243) q[2];
sx q[2];
rz(-2.2218349) q[2];
sx q[2];
rz(0.541614) q[2];
rz(1.2540865) q[3];
sx q[3];
rz(-1.5897635) q[3];
sx q[3];
rz(0.90887535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33614531) q[0];
sx q[0];
rz(-0.31444028) q[0];
sx q[0];
rz(-1.0910777) q[0];
rz(-2.3167141) q[1];
sx q[1];
rz(-0.42366091) q[1];
sx q[1];
rz(-1.0728015) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0925508) q[0];
sx q[0];
rz(-1.9922087) q[0];
sx q[0];
rz(1.7803935) q[0];
rz(1.593441) q[2];
sx q[2];
rz(-1.5839108) q[2];
sx q[2];
rz(0.85044059) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3531216) q[1];
sx q[1];
rz(-1.7040729) q[1];
sx q[1];
rz(-0.2443831) q[1];
rz(-pi) q[2];
rz(1.1680702) q[3];
sx q[3];
rz(-1.793705) q[3];
sx q[3];
rz(0.18833489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23952809) q[2];
sx q[2];
rz(-2.5967279) q[2];
sx q[2];
rz(-5/(6*pi)) q[2];
rz(-0.15051633) q[3];
sx q[3];
rz(-1.1082606) q[3];
sx q[3];
rz(-0.36293852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4758258) q[0];
sx q[0];
rz(-3.0463687) q[0];
sx q[0];
rz(1.6640523) q[0];
rz(0.80329576) q[1];
sx q[1];
rz(-1.3731615) q[1];
sx q[1];
rz(1.0172179) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8724339) q[0];
sx q[0];
rz(-2.9000421) q[0];
sx q[0];
rz(1.0046602) q[0];
rz(-pi) q[1];
x q[1];
rz(2.656865) q[2];
sx q[2];
rz(-0.51186168) q[2];
sx q[2];
rz(0.32527015) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.89722952) q[1];
sx q[1];
rz(-1.23037) q[1];
sx q[1];
rz(-1.8805518) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32633304) q[3];
sx q[3];
rz(-1.2830955) q[3];
sx q[3];
rz(0.74244754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.21641714) q[2];
sx q[2];
rz(-1.567013) q[2];
sx q[2];
rz(-2.4851921) q[2];
rz(2.9260855) q[3];
sx q[3];
rz(-1.4114722) q[3];
sx q[3];
rz(-1.4490674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0553174) q[0];
sx q[0];
rz(-1.3441939) q[0];
sx q[0];
rz(2.2917746) q[0];
rz(2.1814116) q[1];
sx q[1];
rz(-0.4147059) q[1];
sx q[1];
rz(-2.2023831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67341833) q[0];
sx q[0];
rz(-0.8020173) q[0];
sx q[0];
rz(1.9413663) q[0];
rz(-pi) q[1];
rz(-1.5884253) q[2];
sx q[2];
rz(-2.7269502) q[2];
sx q[2];
rz(-1.3585684) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43329217) q[1];
sx q[1];
rz(-0.95656779) q[1];
sx q[1];
rz(2.0310165) q[1];
rz(-2.2978503) q[3];
sx q[3];
rz(-2.4795723) q[3];
sx q[3];
rz(-0.5211422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5592929) q[2];
sx q[2];
rz(-1.7377995) q[2];
sx q[2];
rz(-0.55029184) q[2];
rz(0.127921) q[3];
sx q[3];
rz(-0.13851276) q[3];
sx q[3];
rz(-2.1730455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6561103) q[0];
sx q[0];
rz(-0.67018296) q[0];
sx q[0];
rz(-0.67076587) q[0];
rz(0.35519629) q[1];
sx q[1];
rz(-1.6658446) q[1];
sx q[1];
rz(2.7459941) q[1];
rz(-0.14329362) q[2];
sx q[2];
rz(-1.431965) q[2];
sx q[2];
rz(1.9964249) q[2];
rz(-2.0215423) q[3];
sx q[3];
rz(-1.7786296) q[3];
sx q[3];
rz(2.8854388) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
