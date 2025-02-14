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
rz(-2.812204) q[0];
sx q[0];
rz(-0.63151276) q[0];
sx q[0];
rz(3.1407177) q[0];
rz(-0.64088351) q[1];
sx q[1];
rz(5.2823245) q[1];
sx q[1];
rz(9.7699788) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.810173) q[0];
sx q[0];
rz(-1.6616482) q[0];
sx q[0];
rz(-3.1167555) q[0];
x q[1];
rz(2.7617747) q[2];
sx q[2];
rz(-1.7079412) q[2];
sx q[2];
rz(-1.6452546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.36269128) q[1];
sx q[1];
rz(-0.76299113) q[1];
sx q[1];
rz(0.94339006) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3561068) q[3];
sx q[3];
rz(-0.076138894) q[3];
sx q[3];
rz(-2.3977181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.22902809) q[2];
sx q[2];
rz(-1.3338858) q[2];
sx q[2];
rz(-0.20797569) q[2];
rz(-0.29911706) q[3];
sx q[3];
rz(-0.57204539) q[3];
sx q[3];
rz(-2.0007029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3308554) q[0];
sx q[0];
rz(-2.9006697) q[0];
sx q[0];
rz(-2.148707) q[0];
rz(-1.8244686) q[1];
sx q[1];
rz(-2.8254852) q[1];
sx q[1];
rz(0.77450007) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68607722) q[0];
sx q[0];
rz(-2.9630399) q[0];
sx q[0];
rz(1.5088085) q[0];
x q[1];
rz(-2.3341137) q[2];
sx q[2];
rz(-1.9614855) q[2];
sx q[2];
rz(-1.9964068) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.61770536) q[1];
sx q[1];
rz(-1.7951709) q[1];
sx q[1];
rz(0.38274204) q[1];
x q[2];
rz(-1.1011613) q[3];
sx q[3];
rz(-1.3691994) q[3];
sx q[3];
rz(2.2045709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.896686) q[2];
sx q[2];
rz(-1.2865571) q[2];
sx q[2];
rz(2.1265105) q[2];
rz(-0.24909881) q[3];
sx q[3];
rz(-2.2738012) q[3];
sx q[3];
rz(0.36756137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2772813) q[0];
sx q[0];
rz(-0.69121498) q[0];
sx q[0];
rz(-0.15750289) q[0];
rz(-1.0109673) q[1];
sx q[1];
rz(-2.8087661) q[1];
sx q[1];
rz(-0.13883042) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75051266) q[0];
sx q[0];
rz(-1.8939928) q[0];
sx q[0];
rz(-1.9531519) q[0];
x q[1];
rz(-2.6723318) q[2];
sx q[2];
rz(-0.77232328) q[2];
sx q[2];
rz(-2.0283716) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.638354) q[1];
sx q[1];
rz(-1.8625096) q[1];
sx q[1];
rz(0.20511638) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7458472) q[3];
sx q[3];
rz(-1.850046) q[3];
sx q[3];
rz(-1.4917013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46131721) q[2];
sx q[2];
rz(-2.2291144) q[2];
sx q[2];
rz(-2.7834564) q[2];
rz(-2.4628468) q[3];
sx q[3];
rz(-2.1923784) q[3];
sx q[3];
rz(-1.0391191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045227483) q[0];
sx q[0];
rz(-2.1684833) q[0];
sx q[0];
rz(0.3048234) q[0];
rz(2.7335956) q[1];
sx q[1];
rz(-1.7034737) q[1];
sx q[1];
rz(0.92794424) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.613425) q[0];
sx q[0];
rz(-1.6990464) q[0];
sx q[0];
rz(1.3750465) q[0];
rz(2.1963652) q[2];
sx q[2];
rz(-2.0778227) q[2];
sx q[2];
rz(2.8692109) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9548861) q[1];
sx q[1];
rz(-2.2518603) q[1];
sx q[1];
rz(-1.7658556) q[1];
rz(-pi) q[2];
rz(0.66189142) q[3];
sx q[3];
rz(-2.3958979) q[3];
sx q[3];
rz(-2.6357366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1912332) q[2];
sx q[2];
rz(-0.57098907) q[2];
sx q[2];
rz(-0.18816571) q[2];
rz(-1.865271) q[3];
sx q[3];
rz(-1.8264344) q[3];
sx q[3];
rz(1.7108542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6998049) q[0];
sx q[0];
rz(-2.7975174) q[0];
sx q[0];
rz(3.1222043) q[0];
rz(2.0806606) q[1];
sx q[1];
rz(-1.414199) q[1];
sx q[1];
rz(1.9020938) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083154924) q[0];
sx q[0];
rz(-0.64955901) q[0];
sx q[0];
rz(-2.044528) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4475736) q[2];
sx q[2];
rz(-2.2713714) q[2];
sx q[2];
rz(-1.4365591) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.46249786) q[1];
sx q[1];
rz(-1.0626918) q[1];
sx q[1];
rz(-1.2485882) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.40941671) q[3];
sx q[3];
rz(-1.1077266) q[3];
sx q[3];
rz(0.45209979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0197319) q[2];
sx q[2];
rz(-1.8283565) q[2];
sx q[2];
rz(1.8149553) q[2];
rz(2.895368) q[3];
sx q[3];
rz(-1.1028057) q[3];
sx q[3];
rz(-0.8518014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.6061848) q[0];
sx q[0];
rz(-0.98537213) q[0];
sx q[0];
rz(2.4024409) q[0];
rz(2.7958561) q[1];
sx q[1];
rz(-1.812499) q[1];
sx q[1];
rz(-0.87127042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7428674) q[0];
sx q[0];
rz(-2.350432) q[0];
sx q[0];
rz(1.6683116) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1484457) q[2];
sx q[2];
rz(-0.9075853) q[2];
sx q[2];
rz(0.011485966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.516219) q[1];
sx q[1];
rz(-1.1828848) q[1];
sx q[1];
rz(0.078887786) q[1];
rz(1.096181) q[3];
sx q[3];
rz(-2.1283009) q[3];
sx q[3];
rz(1.4589918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8295916) q[2];
sx q[2];
rz(-2.5025554) q[2];
sx q[2];
rz(-0.87257067) q[2];
rz(1.5729337) q[3];
sx q[3];
rz(-2.5376153) q[3];
sx q[3];
rz(-0.93305552) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0503814) q[0];
sx q[0];
rz(-2.8822883) q[0];
sx q[0];
rz(0.76164371) q[0];
rz(-2.7161982) q[1];
sx q[1];
rz(-2.0014706) q[1];
sx q[1];
rz(2.2147307) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2180397) q[0];
sx q[0];
rz(-1.3430376) q[0];
sx q[0];
rz(-2.8995598) q[0];
rz(-pi) q[1];
rz(2.2144187) q[2];
sx q[2];
rz(-2.1774051) q[2];
sx q[2];
rz(1.7546897) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.35173479) q[1];
sx q[1];
rz(-0.42966336) q[1];
sx q[1];
rz(-1.0142782) q[1];
rz(-0.94759946) q[3];
sx q[3];
rz(-1.5316846) q[3];
sx q[3];
rz(2.2653803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0285792) q[2];
sx q[2];
rz(-0.69602746) q[2];
sx q[2];
rz(-0.87116233) q[2];
rz(2.8431559) q[3];
sx q[3];
rz(-1.2454183) q[3];
sx q[3];
rz(1.8106073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4153862) q[0];
sx q[0];
rz(-0.50006777) q[0];
sx q[0];
rz(2.723208) q[0];
rz(-2.8437974) q[1];
sx q[1];
rz(-1.9458385) q[1];
sx q[1];
rz(-0.13430886) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1100562) q[0];
sx q[0];
rz(-2.2087065) q[0];
sx q[0];
rz(-2.4619034) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18608002) q[2];
sx q[2];
rz(-1.1717516) q[2];
sx q[2];
rz(3.0732791) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.011439104) q[1];
sx q[1];
rz(-1.7053889) q[1];
sx q[1];
rz(2.9573739) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0487324) q[3];
sx q[3];
rz(-0.52937859) q[3];
sx q[3];
rz(0.10806187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40992752) q[2];
sx q[2];
rz(-1.2890559) q[2];
sx q[2];
rz(2.4995787) q[2];
rz(-2.0783453) q[3];
sx q[3];
rz(-0.65170538) q[3];
sx q[3];
rz(-2.1003907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6006271) q[0];
sx q[0];
rz(-0.0890812) q[0];
sx q[0];
rz(-0.11216057) q[0];
rz(-1.3345831) q[1];
sx q[1];
rz(-2.5490675) q[1];
sx q[1];
rz(2.1122011) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0784796) q[0];
sx q[0];
rz(-1.6029583) q[0];
sx q[0];
rz(-0.86634791) q[0];
x q[1];
rz(0.30390988) q[2];
sx q[2];
rz(-1.5329156) q[2];
sx q[2];
rz(1.607995) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9255213) q[1];
sx q[1];
rz(-2.9072484) q[1];
sx q[1];
rz(-1.7611124) q[1];
rz(2.9618046) q[3];
sx q[3];
rz(-0.4745634) q[3];
sx q[3];
rz(-3.0347412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.42334291) q[2];
sx q[2];
rz(-3.0666879) q[2];
sx q[2];
rz(-3.1035799) q[2];
rz(0.612261) q[3];
sx q[3];
rz(-0.93658787) q[3];
sx q[3];
rz(-3.0003701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8925979) q[0];
sx q[0];
rz(-1.4709512) q[0];
sx q[0];
rz(2.7227962) q[0];
rz(1.8839802) q[1];
sx q[1];
rz(-1.5092756) q[1];
sx q[1];
rz(2.9990101) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0973546) q[0];
sx q[0];
rz(-1.729106) q[0];
sx q[0];
rz(-1.5990785) q[0];
rz(-pi) q[1];
rz(-2.5129287) q[2];
sx q[2];
rz(-1.3445879) q[2];
sx q[2];
rz(-2.5369801) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0847454) q[1];
sx q[1];
rz(-2.1386302) q[1];
sx q[1];
rz(2.826087) q[1];
x q[2];
rz(-0.94513388) q[3];
sx q[3];
rz(-0.54900733) q[3];
sx q[3];
rz(1.6459203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3451781) q[2];
sx q[2];
rz(-0.08064457) q[2];
sx q[2];
rz(-0.26819116) q[2];
rz(1.4043407) q[3];
sx q[3];
rz(-2.0495448) q[3];
sx q[3];
rz(-1.0442737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0161229) q[0];
sx q[0];
rz(-2.7653427) q[0];
sx q[0];
rz(0.34633037) q[0];
rz(3.012433) q[1];
sx q[1];
rz(-1.8864514) q[1];
sx q[1];
rz(1.4056978) q[1];
rz(-1.9228946) q[2];
sx q[2];
rz(-1.1209956) q[2];
sx q[2];
rz(1.6605177) q[2];
rz(1.3553452) q[3];
sx q[3];
rz(-0.30134311) q[3];
sx q[3];
rz(-2.338196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
