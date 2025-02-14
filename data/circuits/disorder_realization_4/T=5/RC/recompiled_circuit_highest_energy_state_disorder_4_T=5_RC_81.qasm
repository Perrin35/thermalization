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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1064426) q[0];
sx q[0];
rz(-1.5272836) q[0];
sx q[0];
rz(-1.6405126) q[0];
rz(0.040272399) q[2];
sx q[2];
rz(-1.2408409) q[2];
sx q[2];
rz(1.9391935) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1920584) q[1];
sx q[1];
rz(-0.84334521) q[1];
sx q[1];
rz(0.53967584) q[1];
rz(-pi) q[2];
rz(0.51835097) q[3];
sx q[3];
rz(-0.96203321) q[3];
sx q[3];
rz(1.4639554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4618571) q[2];
sx q[2];
rz(-1.4718141) q[2];
sx q[2];
rz(2.8409345) q[2];
rz(-2.4833208) q[3];
sx q[3];
rz(-2.7418147) q[3];
sx q[3];
rz(-0.58832735) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7268426) q[0];
sx q[0];
rz(-0.86597935) q[0];
sx q[0];
rz(-0.17587371) q[0];
rz(2.580592) q[1];
sx q[1];
rz(-2.439552) q[1];
sx q[1];
rz(-2.9452513) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25717673) q[0];
sx q[0];
rz(-1.1286583) q[0];
sx q[0];
rz(2.7701562) q[0];
rz(-pi) q[1];
rz(-3.0220993) q[2];
sx q[2];
rz(-1.4548848) q[2];
sx q[2];
rz(-1.8279861) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2113116) q[1];
sx q[1];
rz(-2.086211) q[1];
sx q[1];
rz(-3.1392908) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34336932) q[3];
sx q[3];
rz(-1.6968622) q[3];
sx q[3];
rz(2.8186563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.089513) q[2];
sx q[2];
rz(-2.5255346) q[2];
sx q[2];
rz(1.4698131) q[2];
rz(0.52516627) q[3];
sx q[3];
rz(-1.4273806) q[3];
sx q[3];
rz(-0.66450459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1216275) q[0];
sx q[0];
rz(-2.2522734) q[0];
sx q[0];
rz(0.63051939) q[0];
rz(1.1446674) q[1];
sx q[1];
rz(-1.8960709) q[1];
sx q[1];
rz(3.0847881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5543514) q[0];
sx q[0];
rz(-1.2040753) q[0];
sx q[0];
rz(1.1318867) q[0];
rz(-pi) q[1];
rz(0.47614758) q[2];
sx q[2];
rz(-0.94149977) q[2];
sx q[2];
rz(-2.4193633) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5377429) q[1];
sx q[1];
rz(-1.7122388) q[1];
sx q[1];
rz(-0.46700041) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1435236) q[3];
sx q[3];
rz(-1.3766854) q[3];
sx q[3];
rz(-2.8097092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28222617) q[2];
sx q[2];
rz(-1.4713918) q[2];
sx q[2];
rz(2.7024506) q[2];
rz(-1.7450843) q[3];
sx q[3];
rz(-1.8789995) q[3];
sx q[3];
rz(2.5631574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0466995) q[0];
sx q[0];
rz(-0.72143227) q[0];
sx q[0];
rz(-2.1266345) q[0];
rz(-0.062648423) q[1];
sx q[1];
rz(-0.95476127) q[1];
sx q[1];
rz(0.28422022) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009722) q[0];
sx q[0];
rz(-0.56864244) q[0];
sx q[0];
rz(2.4533662) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7106871) q[2];
sx q[2];
rz(-2.0645185) q[2];
sx q[2];
rz(1.7715275) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41955796) q[1];
sx q[1];
rz(-2.7830208) q[1];
sx q[1];
rz(0.034073985) q[1];
rz(-1.5169859) q[3];
sx q[3];
rz(-2.887886) q[3];
sx q[3];
rz(-2.1068609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8404428) q[2];
sx q[2];
rz(-2.132685) q[2];
sx q[2];
rz(-2.3811049) q[2];
rz(2.8216951) q[3];
sx q[3];
rz(-1.1287929) q[3];
sx q[3];
rz(1.0285876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496721) q[0];
sx q[0];
rz(-0.55067647) q[0];
sx q[0];
rz(0.40529761) q[0];
rz(2.6536476) q[1];
sx q[1];
rz(-1.1064203) q[1];
sx q[1];
rz(2.778756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3317446) q[0];
sx q[0];
rz(-2.6001626) q[0];
sx q[0];
rz(0.62744139) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9486729) q[2];
sx q[2];
rz(-1.0185756) q[2];
sx q[2];
rz(0.044980031) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0367485) q[1];
sx q[1];
rz(-1.3499773) q[1];
sx q[1];
rz(-1.6206431) q[1];
rz(-0.36146253) q[3];
sx q[3];
rz(-0.72849724) q[3];
sx q[3];
rz(-2.0036774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.001658) q[2];
sx q[2];
rz(-1.1589061) q[2];
sx q[2];
rz(0.72251594) q[2];
rz(-2.6288988) q[3];
sx q[3];
rz(-2.2721458) q[3];
sx q[3];
rz(-1.2041913) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743415) q[0];
sx q[0];
rz(-2.5047472) q[0];
sx q[0];
rz(-0.27477086) q[0];
rz(-0.35708669) q[1];
sx q[1];
rz(-1.5029224) q[1];
sx q[1];
rz(-1.9836609) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9102218) q[0];
sx q[0];
rz(-1.3274258) q[0];
sx q[0];
rz(-2.2621808) q[0];
rz(-0.91041301) q[2];
sx q[2];
rz(-2.0636301) q[2];
sx q[2];
rz(0.016066859) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14640644) q[1];
sx q[1];
rz(-0.9764834) q[1];
sx q[1];
rz(1.8176778) q[1];
rz(-pi) q[2];
rz(2.4729608) q[3];
sx q[3];
rz(-0.9210081) q[3];
sx q[3];
rz(-1.970552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5623515) q[2];
sx q[2];
rz(-1.3902105) q[2];
sx q[2];
rz(-2.5518899) q[2];
rz(-1.3127182) q[3];
sx q[3];
rz(-1.4785942) q[3];
sx q[3];
rz(2.5780799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8967459) q[0];
sx q[0];
rz(-0.99269301) q[0];
sx q[0];
rz(2.8712811) q[0];
rz(2.7145794) q[1];
sx q[1];
rz(-1.7662363) q[1];
sx q[1];
rz(0.40831533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9664551) q[0];
sx q[0];
rz(-1.2027367) q[0];
sx q[0];
rz(1.4428663) q[0];
rz(-3.0403942) q[2];
sx q[2];
rz(-2.2658341) q[2];
sx q[2];
rz(-1.1388495) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.282849) q[1];
sx q[1];
rz(-0.74902422) q[1];
sx q[1];
rz(0.94094679) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72260489) q[3];
sx q[3];
rz(-1.4935023) q[3];
sx q[3];
rz(2.3376436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.938505) q[2];
sx q[2];
rz(-1.0838584) q[2];
sx q[2];
rz(-0.92424029) q[2];
rz(0.83485323) q[3];
sx q[3];
rz(-0.82930851) q[3];
sx q[3];
rz(0.99669641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95558178) q[0];
sx q[0];
rz(-2.7614433) q[0];
sx q[0];
rz(-0.36744776) q[0];
rz(-0.5683178) q[1];
sx q[1];
rz(-1.808993) q[1];
sx q[1];
rz(2.7265991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58632219) q[0];
sx q[0];
rz(-1.8316226) q[0];
sx q[0];
rz(2.0404979) q[0];
rz(0.80525406) q[2];
sx q[2];
rz(-0.70119154) q[2];
sx q[2];
rz(2.5561577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.69232443) q[1];
sx q[1];
rz(-1.0639166) q[1];
sx q[1];
rz(-2.436422) q[1];
rz(-pi) q[2];
rz(-1.5852195) q[3];
sx q[3];
rz(-0.34581772) q[3];
sx q[3];
rz(-2.7399969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7830398) q[2];
sx q[2];
rz(-0.8873322) q[2];
sx q[2];
rz(-2.3084194) q[2];
rz(2.7821275) q[3];
sx q[3];
rz(-2.0514252) q[3];
sx q[3];
rz(1.5045213) q[3];
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
rz(-pi/2) q[3];
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
rz(1.1737162) q[0];
sx q[0];
rz(-2.739527) q[0];
sx q[0];
rz(-0.014884431) q[0];
rz(1.124565) q[1];
sx q[1];
rz(-2.896307) q[1];
sx q[1];
rz(3.0344149) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8421321) q[0];
sx q[0];
rz(-0.74343527) q[0];
sx q[0];
rz(-1.7920262) q[0];
x q[1];
rz(0.46469073) q[2];
sx q[2];
rz(-2.26311) q[2];
sx q[2];
rz(-1.2274881) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8228777) q[1];
sx q[1];
rz(-1.4221622) q[1];
sx q[1];
rz(-2.6591402) q[1];
rz(0.26132432) q[3];
sx q[3];
rz(-1.7423986) q[3];
sx q[3];
rz(1.8971407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5643481) q[2];
sx q[2];
rz(-1.8337092) q[2];
sx q[2];
rz(2.2819819) q[2];
rz(2.882242) q[3];
sx q[3];
rz(-1.7813464) q[3];
sx q[3];
rz(0.41659659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49941007) q[0];
sx q[0];
rz(-0.12633093) q[0];
sx q[0];
rz(-1.1580178) q[0];
rz(0.49531373) q[1];
sx q[1];
rz(-1.1664349) q[1];
sx q[1];
rz(-2.8445629) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3619694) q[0];
sx q[0];
rz(-1.2517778) q[0];
sx q[0];
rz(-3.1373128) q[0];
rz(3.0115602) q[2];
sx q[2];
rz(-2.8388925) q[2];
sx q[2];
rz(-1.5695733) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.90290711) q[1];
sx q[1];
rz(-2.8082962) q[1];
sx q[1];
rz(-0.500228) q[1];
x q[2];
rz(-0.90822753) q[3];
sx q[3];
rz(-1.5836827) q[3];
sx q[3];
rz(2.9290309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.872252) q[2];
sx q[2];
rz(-0.38186914) q[2];
sx q[2];
rz(-2.2807109) q[2];
rz(-2.6779029) q[3];
sx q[3];
rz(-2.2380565) q[3];
sx q[3];
rz(-2.004682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1562445) q[0];
sx q[0];
rz(-1.5736268) q[0];
sx q[0];
rz(-1.9851984) q[0];
rz(1.3337749) q[1];
sx q[1];
rz(-1.6009686) q[1];
sx q[1];
rz(0.70645465) q[1];
rz(2.4598224) q[2];
sx q[2];
rz(-0.78777704) q[2];
sx q[2];
rz(-2.4252979) q[2];
rz(1.9577338) q[3];
sx q[3];
rz(-2.0328184) q[3];
sx q[3];
rz(0.55749374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
