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
rz(-1.8972549) q[0];
sx q[0];
rz(-2.3490348) q[0];
sx q[0];
rz(-2.8886524) q[0];
rz(5.5090299) q[1];
sx q[1];
rz(2.7634662) q[1];
sx q[1];
rz(9.0378349) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88322631) q[0];
sx q[0];
rz(-2.7233426) q[0];
sx q[0];
rz(0.89624385) q[0];
x q[1];
rz(-2.5637881) q[2];
sx q[2];
rz(-1.6434323) q[2];
sx q[2];
rz(-2.5355123) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7751037) q[1];
sx q[1];
rz(-1.4867884) q[1];
sx q[1];
rz(-1.4909527) q[1];
x q[2];
rz(2.0617332) q[3];
sx q[3];
rz(-2.1055974) q[3];
sx q[3];
rz(0.25568572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0953377) q[2];
sx q[2];
rz(-2.3825808) q[2];
sx q[2];
rz(-1.7389899) q[2];
rz(1.4143573) q[3];
sx q[3];
rz(-0.71310133) q[3];
sx q[3];
rz(-0.49899092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801179) q[0];
sx q[0];
rz(-2.6533227) q[0];
sx q[0];
rz(-0.71037355) q[0];
rz(-1.0164227) q[1];
sx q[1];
rz(-1.8417532) q[1];
sx q[1];
rz(-2.3764835) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.13694) q[0];
sx q[0];
rz(-1.0400794) q[0];
sx q[0];
rz(-1.9416481) q[0];
x q[1];
rz(0.76025195) q[2];
sx q[2];
rz(-2.2905802) q[2];
sx q[2];
rz(-2.8362897) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.930508) q[1];
sx q[1];
rz(-0.74568664) q[1];
sx q[1];
rz(2.0603643) q[1];
x q[2];
rz(-1.2364976) q[3];
sx q[3];
rz(-1.7124885) q[3];
sx q[3];
rz(0.96135512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1171099) q[2];
sx q[2];
rz(-2.4958002) q[2];
sx q[2];
rz(-2.4369241) q[2];
rz(2.9674528) q[3];
sx q[3];
rz(-0.87248674) q[3];
sx q[3];
rz(2.7599938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0524549) q[0];
sx q[0];
rz(-1.8280886) q[0];
sx q[0];
rz(-2.254159) q[0];
rz(-1.8916091) q[1];
sx q[1];
rz(-1.9307815) q[1];
sx q[1];
rz(-1.3608305) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75653532) q[0];
sx q[0];
rz(-0.65246118) q[0];
sx q[0];
rz(0.65076179) q[0];
rz(-pi) q[1];
rz(3.0615882) q[2];
sx q[2];
rz(-1.4334049) q[2];
sx q[2];
rz(-0.72890711) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8831142) q[1];
sx q[1];
rz(-0.7683903) q[1];
sx q[1];
rz(1.7151296) q[1];
rz(-pi) q[2];
rz(-1.181608) q[3];
sx q[3];
rz(-1.5312315) q[3];
sx q[3];
rz(-0.54007441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7736194) q[2];
sx q[2];
rz(-0.36586389) q[2];
sx q[2];
rz(1.492929) q[2];
rz(-0.44935539) q[3];
sx q[3];
rz(-1.2949233) q[3];
sx q[3];
rz(-1.2202643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0538977) q[0];
sx q[0];
rz(-0.59006417) q[0];
sx q[0];
rz(1.4087403) q[0];
rz(0.98211163) q[1];
sx q[1];
rz(-1.5928007) q[1];
sx q[1];
rz(-1.7353479) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8187253) q[0];
sx q[0];
rz(-1.3565085) q[0];
sx q[0];
rz(0.11361285) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1020145) q[2];
sx q[2];
rz(-1.5180032) q[2];
sx q[2];
rz(0.18597183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50735578) q[1];
sx q[1];
rz(-1.0546396) q[1];
sx q[1];
rz(0.28395758) q[1];
x q[2];
rz(2.2194524) q[3];
sx q[3];
rz(-1.7867844) q[3];
sx q[3];
rz(0.43555015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1534319) q[2];
sx q[2];
rz(-2.1250696) q[2];
sx q[2];
rz(-2.9856258) q[2];
rz(-2.4912452) q[3];
sx q[3];
rz(-1.1740843) q[3];
sx q[3];
rz(0.11421886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0735737) q[0];
sx q[0];
rz(-1.2696215) q[0];
sx q[0];
rz(-2.9408348) q[0];
rz(-1.1772032) q[1];
sx q[1];
rz(-2.2889844) q[1];
sx q[1];
rz(-1.1600561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5302488) q[0];
sx q[0];
rz(-0.28168105) q[0];
sx q[0];
rz(0.92667093) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.54520901) q[2];
sx q[2];
rz(-0.13449796) q[2];
sx q[2];
rz(-0.70869499) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9721891) q[1];
sx q[1];
rz(-2.6693925) q[1];
sx q[1];
rz(2.4257823) q[1];
x q[2];
rz(1.2358642) q[3];
sx q[3];
rz(-2.5865002) q[3];
sx q[3];
rz(0.60716437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3843627) q[2];
sx q[2];
rz(-2.195916) q[2];
sx q[2];
rz(-0.45822701) q[2];
rz(2.5107757) q[3];
sx q[3];
rz(-1.2354555) q[3];
sx q[3];
rz(0.49921504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.683627) q[0];
sx q[0];
rz(-2.2891335) q[0];
sx q[0];
rz(3.1347347) q[0];
rz(-0.231617) q[1];
sx q[1];
rz(-1.3791142) q[1];
sx q[1];
rz(-1.0622271) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4214913) q[0];
sx q[0];
rz(-1.4473697) q[0];
sx q[0];
rz(-2.8272259) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49906667) q[2];
sx q[2];
rz(-1.7573414) q[2];
sx q[2];
rz(1.3908433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66328996) q[1];
sx q[1];
rz(-1.3172825) q[1];
sx q[1];
rz(-0.67914825) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3788578) q[3];
sx q[3];
rz(-0.71832359) q[3];
sx q[3];
rz(3.0305221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0007533) q[2];
sx q[2];
rz(-1.8893628) q[2];
sx q[2];
rz(1.2678649) q[2];
rz(-2.2800692) q[3];
sx q[3];
rz(-1.0637161) q[3];
sx q[3];
rz(-1.4388194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.2086901) q[0];
sx q[0];
rz(-0.33354315) q[0];
sx q[0];
rz(2.9260337) q[0];
rz(0.49628273) q[1];
sx q[1];
rz(-1.9535306) q[1];
sx q[1];
rz(0.15942474) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8474583) q[0];
sx q[0];
rz(-2.205664) q[0];
sx q[0];
rz(-3.1383187) q[0];
rz(-0.36783571) q[2];
sx q[2];
rz(-2.6421037) q[2];
sx q[2];
rz(-2.6373581) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9558952) q[1];
sx q[1];
rz(-1.5296827) q[1];
sx q[1];
rz(-1.3354882) q[1];
rz(-1.0219021) q[3];
sx q[3];
rz(-1.1682142) q[3];
sx q[3];
rz(2.2233913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40519199) q[2];
sx q[2];
rz(-2.0173732) q[2];
sx q[2];
rz(-2.6250725) q[2];
rz(1.9308331) q[3];
sx q[3];
rz(-1.4529994) q[3];
sx q[3];
rz(1.7621382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(1.6571758) q[0];
sx q[0];
rz(-0.076198904) q[0];
sx q[0];
rz(-2.6013689) q[0];
rz(-1.5243439) q[1];
sx q[1];
rz(-1.0018307) q[1];
sx q[1];
rz(1.8744972) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0277242) q[0];
sx q[0];
rz(-1.0719711) q[0];
sx q[0];
rz(-1.3214825) q[0];
rz(0.71155602) q[2];
sx q[2];
rz(-0.91469736) q[2];
sx q[2];
rz(3.1392821) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.46173501) q[1];
sx q[1];
rz(-1.2038132) q[1];
sx q[1];
rz(1.4162029) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2121939) q[3];
sx q[3];
rz(-1.2817973) q[3];
sx q[3];
rz(-2.3232587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2544864) q[2];
sx q[2];
rz(-1.9820513) q[2];
sx q[2];
rz(2.9642588) q[2];
rz(-1.0698498) q[3];
sx q[3];
rz(-2.0900574) q[3];
sx q[3];
rz(0.18226084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3607445) q[0];
sx q[0];
rz(-1.9875263) q[0];
sx q[0];
rz(0.10072197) q[0];
rz(1.992647) q[1];
sx q[1];
rz(-1.9457685) q[1];
sx q[1];
rz(-1.1740059) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74885313) q[0];
sx q[0];
rz(-2.1707188) q[0];
sx q[0];
rz(-1.8003182) q[0];
rz(-0.76624932) q[2];
sx q[2];
rz(-1.6403997) q[2];
sx q[2];
rz(2.0947667) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29193297) q[1];
sx q[1];
rz(-1.5330452) q[1];
sx q[1];
rz(-3.1221366) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2618746) q[3];
sx q[3];
rz(-1.4914601) q[3];
sx q[3];
rz(0.21896958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.164087) q[2];
sx q[2];
rz(-1.6567433) q[2];
sx q[2];
rz(-0.40317765) q[2];
rz(0.66344231) q[3];
sx q[3];
rz(-0.57938975) q[3];
sx q[3];
rz(-0.0042393953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73096257) q[0];
sx q[0];
rz(-1.0799438) q[0];
sx q[0];
rz(0.078911111) q[0];
rz(-2.2321189) q[1];
sx q[1];
rz(-2.0957004) q[1];
sx q[1];
rz(-0.39631072) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3672597) q[0];
sx q[0];
rz(-3.1235187) q[0];
sx q[0];
rz(2.0967183) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0070856241) q[2];
sx q[2];
rz(-2.3478697) q[2];
sx q[2];
rz(-0.34863472) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2036613) q[1];
sx q[1];
rz(-1.0962558) q[1];
sx q[1];
rz(-1.306626) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86045281) q[3];
sx q[3];
rz(-2.3858968) q[3];
sx q[3];
rz(-2.6967654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.69184715) q[2];
sx q[2];
rz(-2.2578466) q[2];
sx q[2];
rz(0.79908243) q[2];
rz(-3.0633022) q[3];
sx q[3];
rz(-1.8872063) q[3];
sx q[3];
rz(-2.4962795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66452022) q[0];
sx q[0];
rz(-1.399853) q[0];
sx q[0];
rz(0.54367263) q[0];
rz(1.1935344) q[1];
sx q[1];
rz(-1.8652893) q[1];
sx q[1];
rz(0.51020772) q[1];
rz(1.1134182) q[2];
sx q[2];
rz(-1.7241679) q[2];
sx q[2];
rz(0.67571251) q[2];
rz(2.5922248) q[3];
sx q[3];
rz(-2.3372169) q[3];
sx q[3];
rz(-1.611471) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
