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
rz(0.34613553) q[0];
sx q[0];
rz(4.7799568) q[0];
sx q[0];
rz(10.156375) q[0];
rz(-5.8335891) q[1];
sx q[1];
rz(3.3291266) q[1];
sx q[1];
rz(13.000288) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2303378) q[0];
sx q[0];
rz(-0.83090913) q[0];
sx q[0];
rz(-0.047538443) q[0];
rz(0.97831867) q[2];
sx q[2];
rz(-1.0772395) q[2];
sx q[2];
rz(1.6101642) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.190769) q[1];
sx q[1];
rz(-2.1956823) q[1];
sx q[1];
rz(-0.48887078) q[1];
x q[2];
rz(1.3796666) q[3];
sx q[3];
rz(-1.6014631) q[3];
sx q[3];
rz(2.994184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.71243858) q[2];
sx q[2];
rz(-1.466208) q[2];
sx q[2];
rz(1.0038556) q[2];
rz(0.53576523) q[3];
sx q[3];
rz(-0.86447132) q[3];
sx q[3];
rz(0.0011477688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67343229) q[0];
sx q[0];
rz(-0.68244451) q[0];
sx q[0];
rz(-0.72714192) q[0];
rz(3.0739821) q[1];
sx q[1];
rz(-1.7865684) q[1];
sx q[1];
rz(-2.0179857) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7794117) q[0];
sx q[0];
rz(-1.0955278) q[0];
sx q[0];
rz(-1.2527466) q[0];
rz(3.0762663) q[2];
sx q[2];
rz(-1.4484754) q[2];
sx q[2];
rz(2.9768012) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3930291) q[1];
sx q[1];
rz(-1.2911699) q[1];
sx q[1];
rz(0.10348998) q[1];
x q[2];
rz(-2.0455849) q[3];
sx q[3];
rz(-1.1752306) q[3];
sx q[3];
rz(2.6717348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8565389) q[2];
sx q[2];
rz(-1.5679789) q[2];
sx q[2];
rz(-0.48173586) q[2];
rz(-0.5528062) q[3];
sx q[3];
rz(-1.0779287) q[3];
sx q[3];
rz(-3.0520181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8168617) q[0];
sx q[0];
rz(-2.6664) q[0];
sx q[0];
rz(0.09356308) q[0];
rz(-1.7612673) q[1];
sx q[1];
rz(-0.9535791) q[1];
sx q[1];
rz(2.5043452) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9910332) q[0];
sx q[0];
rz(-2.4399418) q[0];
sx q[0];
rz(2.4392564) q[0];
rz(1.797051) q[2];
sx q[2];
rz(-2.166894) q[2];
sx q[2];
rz(1.4486194) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90594572) q[1];
sx q[1];
rz(-1.596758) q[1];
sx q[1];
rz(-1.173255) q[1];
x q[2];
rz(-0.29070694) q[3];
sx q[3];
rz(-1.6120496) q[3];
sx q[3];
rz(-2.8009529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9351585) q[2];
sx q[2];
rz(-0.55525246) q[2];
sx q[2];
rz(-0.68378249) q[2];
rz(0.23009662) q[3];
sx q[3];
rz(-1.4068539) q[3];
sx q[3];
rz(2.7195215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.049659599) q[0];
sx q[0];
rz(-2.6347418) q[0];
sx q[0];
rz(-1.4403213) q[0];
rz(-0.47538844) q[1];
sx q[1];
rz(-2.5071867) q[1];
sx q[1];
rz(2.5604274) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6027045) q[0];
sx q[0];
rz(-1.4853444) q[0];
sx q[0];
rz(-0.84792015) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0175242) q[2];
sx q[2];
rz(-0.5047732) q[2];
sx q[2];
rz(-3.1044132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6254723) q[1];
sx q[1];
rz(-0.56159084) q[1];
sx q[1];
rz(-0.59188868) q[1];
rz(1.7293594) q[3];
sx q[3];
rz(-1.0845636) q[3];
sx q[3];
rz(2.9610046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0541957) q[2];
sx q[2];
rz(-2.9554458) q[2];
sx q[2];
rz(2.6206214) q[2];
rz(-1.4969131) q[3];
sx q[3];
rz(-1.729634) q[3];
sx q[3];
rz(-1.514667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5882551) q[0];
sx q[0];
rz(-0.27565685) q[0];
sx q[0];
rz(-1.1302554) q[0];
rz(-2.0626119) q[1];
sx q[1];
rz(-1.1993473) q[1];
sx q[1];
rz(-1.1111396) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1905907) q[0];
sx q[0];
rz(-1.4050806) q[0];
sx q[0];
rz(2.034305) q[0];
x q[1];
rz(-0.07569282) q[2];
sx q[2];
rz(-0.53097979) q[2];
sx q[2];
rz(-0.39226433) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6592404) q[1];
sx q[1];
rz(-2.2063046) q[1];
sx q[1];
rz(-0.19584943) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54773982) q[3];
sx q[3];
rz(-0.99036694) q[3];
sx q[3];
rz(1.86059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7378716) q[2];
sx q[2];
rz(-0.51509866) q[2];
sx q[2];
rz(-1.0408638) q[2];
rz(0.8199842) q[3];
sx q[3];
rz(-1.2330202) q[3];
sx q[3];
rz(2.5177054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2417004) q[0];
sx q[0];
rz(-2.5427759) q[0];
sx q[0];
rz(-2.4851121) q[0];
rz(1.7680602) q[1];
sx q[1];
rz(-0.7936002) q[1];
sx q[1];
rz(-1.1318644) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0431932) q[0];
sx q[0];
rz(-0.88623673) q[0];
sx q[0];
rz(1.3149904) q[0];
rz(-pi) q[1];
rz(2.1186351) q[2];
sx q[2];
rz(-0.96998397) q[2];
sx q[2];
rz(3.09077) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.060185) q[1];
sx q[1];
rz(-2.2814732) q[1];
sx q[1];
rz(-3.1050472) q[1];
rz(0.15040654) q[3];
sx q[3];
rz(-1.5664706) q[3];
sx q[3];
rz(1.195601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6953096) q[2];
sx q[2];
rz(-2.535847) q[2];
sx q[2];
rz(0.31965762) q[2];
rz(0.22932912) q[3];
sx q[3];
rz(-1.0234443) q[3];
sx q[3];
rz(2.2314609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0730154) q[0];
sx q[0];
rz(-1.1443161) q[0];
sx q[0];
rz(-1.2314433) q[0];
rz(0.07864174) q[1];
sx q[1];
rz(-1.6784724) q[1];
sx q[1];
rz(0.15422779) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.820136) q[0];
sx q[0];
rz(-2.7559782) q[0];
sx q[0];
rz(-0.48954757) q[0];
rz(-1.1604105) q[2];
sx q[2];
rz(-1.2912116) q[2];
sx q[2];
rz(1.9094163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.73772506) q[1];
sx q[1];
rz(-1.939865) q[1];
sx q[1];
rz(-2.6253659) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46040543) q[3];
sx q[3];
rz(-2.6889927) q[3];
sx q[3];
rz(2.3456338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.33588931) q[2];
sx q[2];
rz(-1.5259537) q[2];
sx q[2];
rz(-1.4914782) q[2];
rz(1.7283745) q[3];
sx q[3];
rz(-1.3946984) q[3];
sx q[3];
rz(-1.570805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0989646) q[0];
sx q[0];
rz(-2.0476116) q[0];
sx q[0];
rz(1.6866823) q[0];
rz(-0.97575724) q[1];
sx q[1];
rz(-1.059633) q[1];
sx q[1];
rz(1.3386493) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0473235) q[0];
sx q[0];
rz(-1.9587095) q[0];
sx q[0];
rz(2.3045425) q[0];
rz(-2.0744616) q[2];
sx q[2];
rz(-1.8657473) q[2];
sx q[2];
rz(-0.3857715) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6852831) q[1];
sx q[1];
rz(-1.7830308) q[1];
sx q[1];
rz(-2.2077399) q[1];
x q[2];
rz(0.66338952) q[3];
sx q[3];
rz(-1.8691392) q[3];
sx q[3];
rz(1.364352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51155382) q[2];
sx q[2];
rz(-2.126882) q[2];
sx q[2];
rz(-0.83941984) q[2];
rz(-1.9073585) q[3];
sx q[3];
rz(-2.0525565) q[3];
sx q[3];
rz(-0.031410005) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79427528) q[0];
sx q[0];
rz(-1.3823771) q[0];
sx q[0];
rz(-2.7224702) q[0];
rz(1.6436815) q[1];
sx q[1];
rz(-1.7828015) q[1];
sx q[1];
rz(0.43325123) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.142665) q[0];
sx q[0];
rz(-1.4523399) q[0];
sx q[0];
rz(-2.0303557) q[0];
x q[1];
rz(3.1243043) q[2];
sx q[2];
rz(-2.4845036) q[2];
sx q[2];
rz(-1.2430056) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0440694) q[1];
sx q[1];
rz(-1.8533924) q[1];
sx q[1];
rz(1.0364012) q[1];
rz(-2.562996) q[3];
sx q[3];
rz(-1.3893327) q[3];
sx q[3];
rz(-2.6058448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7822632) q[2];
sx q[2];
rz(-1.1507582) q[2];
sx q[2];
rz(-2.9456054) q[2];
rz(1.1228784) q[3];
sx q[3];
rz(-0.71176088) q[3];
sx q[3];
rz(-1.6897197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48542431) q[0];
sx q[0];
rz(-2.4805785) q[0];
sx q[0];
rz(2.0458903) q[0];
rz(0.51086673) q[1];
sx q[1];
rz(-2.8725862) q[1];
sx q[1];
rz(2.9313472) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4097262) q[0];
sx q[0];
rz(-1.3748226) q[0];
sx q[0];
rz(-0.98742296) q[0];
rz(-3.0265507) q[2];
sx q[2];
rz(-0.91918531) q[2];
sx q[2];
rz(-0.63331214) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89344674) q[1];
sx q[1];
rz(-2.4706744) q[1];
sx q[1];
rz(0.72709658) q[1];
x q[2];
rz(-1.1594335) q[3];
sx q[3];
rz(-2.1307011) q[3];
sx q[3];
rz(-0.066889569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0290587) q[2];
sx q[2];
rz(-1.5052648) q[2];
sx q[2];
rz(-2.0564334) q[2];
rz(-0.30760136) q[3];
sx q[3];
rz(-0.31809536) q[3];
sx q[3];
rz(-0.65348452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45162421) q[0];
sx q[0];
rz(-1.987048) q[0];
sx q[0];
rz(-1.2183627) q[0];
rz(2.4901509) q[1];
sx q[1];
rz(-0.94964288) q[1];
sx q[1];
rz(-2.3655187) q[1];
rz(-3.1262911) q[2];
sx q[2];
rz(-2.2616784) q[2];
sx q[2];
rz(-0.26605284) q[2];
rz(-0.93919803) q[3];
sx q[3];
rz(-1.152335) q[3];
sx q[3];
rz(-1.4622968) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
