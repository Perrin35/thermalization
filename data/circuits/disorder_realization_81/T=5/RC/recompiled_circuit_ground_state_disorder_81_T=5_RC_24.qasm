OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.12299744) q[0];
sx q[0];
rz(-1.7641492) q[0];
sx q[0];
rz(-2.7121845) q[0];
rz(1.6027066) q[1];
sx q[1];
rz(7.9908854) q[1];
sx q[1];
rz(7.8828852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7943952) q[0];
sx q[0];
rz(-1.3566249) q[0];
sx q[0];
rz(1.2783575) q[0];
rz(-pi) q[1];
rz(-0.23878204) q[2];
sx q[2];
rz(-2.9936643) q[2];
sx q[2];
rz(-1.5996358) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0059389) q[1];
sx q[1];
rz(-2.7126813) q[1];
sx q[1];
rz(-3.0687544) q[1];
x q[2];
rz(0.36088636) q[3];
sx q[3];
rz(-1.8461707) q[3];
sx q[3];
rz(2.9831246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8371007) q[2];
sx q[2];
rz(-1.9383177) q[2];
sx q[2];
rz(-2.4486747) q[2];
rz(0.20017008) q[3];
sx q[3];
rz(-0.19014159) q[3];
sx q[3];
rz(2.0076803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.8189341) q[0];
sx q[0];
rz(-0.19911961) q[0];
sx q[0];
rz(1.06485) q[0];
rz(-0.62243593) q[1];
sx q[1];
rz(-1.7935926) q[1];
sx q[1];
rz(2.650824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3773133) q[0];
sx q[0];
rz(-1.547033) q[0];
sx q[0];
rz(1.5615433) q[0];
x q[1];
rz(-0.46892898) q[2];
sx q[2];
rz(-2.6750419) q[2];
sx q[2];
rz(-2.5836437) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35247861) q[1];
sx q[1];
rz(-1.4755469) q[1];
sx q[1];
rz(-0.75525166) q[1];
rz(-1.2964046) q[3];
sx q[3];
rz(-2.9657722) q[3];
sx q[3];
rz(-2.4309513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.48851594) q[2];
sx q[2];
rz(-1.4986897) q[2];
sx q[2];
rz(-0.24031362) q[2];
rz(0.49992418) q[3];
sx q[3];
rz(-2.5559055) q[3];
sx q[3];
rz(-1.2691931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2568473) q[0];
sx q[0];
rz(-0.95504967) q[0];
sx q[0];
rz(-2.1599059) q[0];
rz(0.49199545) q[1];
sx q[1];
rz(-2.056608) q[1];
sx q[1];
rz(-1.5031776) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7273151) q[0];
sx q[0];
rz(-2.2637667) q[0];
sx q[0];
rz(2.5056865) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0122224) q[2];
sx q[2];
rz(-2.7384842) q[2];
sx q[2];
rz(-2.108824) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.43736378) q[1];
sx q[1];
rz(-1.2384602) q[1];
sx q[1];
rz(-1.9802914) q[1];
rz(-pi) q[2];
rz(1.7812626) q[3];
sx q[3];
rz(-0.099848824) q[3];
sx q[3];
rz(-2.9304402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4846399) q[2];
sx q[2];
rz(-2.6159365) q[2];
sx q[2];
rz(-3.0204115) q[2];
rz(-2.1335404) q[3];
sx q[3];
rz(-1.1153778) q[3];
sx q[3];
rz(-1.0813659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0693531) q[0];
sx q[0];
rz(-0.00088748137) q[0];
sx q[0];
rz(0.51373154) q[0];
rz(-0.95798245) q[1];
sx q[1];
rz(-1.999141) q[1];
sx q[1];
rz(1.4428008) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4289249) q[0];
sx q[0];
rz(-2.6972572) q[0];
sx q[0];
rz(-0.39074583) q[0];
rz(-2.6110177) q[2];
sx q[2];
rz(-1.7624035) q[2];
sx q[2];
rz(0.2577739) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5123547) q[1];
sx q[1];
rz(-1.4886453) q[1];
sx q[1];
rz(-2.1869529) q[1];
rz(0.12421457) q[3];
sx q[3];
rz(-2.476177) q[3];
sx q[3];
rz(1.6197325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.55502597) q[2];
sx q[2];
rz(-2.1777966) q[2];
sx q[2];
rz(-1.2653992) q[2];
rz(-1.966656) q[3];
sx q[3];
rz(-0.73052162) q[3];
sx q[3];
rz(-1.9780212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2555399) q[0];
sx q[0];
rz(-2.9386254) q[0];
sx q[0];
rz(-1.8170005) q[0];
rz(-0.96877226) q[1];
sx q[1];
rz(-1.4854393) q[1];
sx q[1];
rz(1.0383777) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2823855) q[0];
sx q[0];
rz(-0.41658724) q[0];
sx q[0];
rz(0.4075012) q[0];
rz(-2.3164301) q[2];
sx q[2];
rz(-2.1405915) q[2];
sx q[2];
rz(-1.2078326) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6250302) q[1];
sx q[1];
rz(-1.5296401) q[1];
sx q[1];
rz(-2.7709318) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65870993) q[3];
sx q[3];
rz(-1.1532056) q[3];
sx q[3];
rz(0.61826784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6394627) q[2];
sx q[2];
rz(-0.30439964) q[2];
sx q[2];
rz(2.2501865) q[2];
rz(1.2616448) q[3];
sx q[3];
rz(-1.5104048) q[3];
sx q[3];
rz(-0.6663028) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0618184) q[0];
sx q[0];
rz(-2.6462055) q[0];
sx q[0];
rz(2.1451982) q[0];
rz(-2.2415316) q[1];
sx q[1];
rz(-0.78949094) q[1];
sx q[1];
rz(-2.9642504) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7686488) q[0];
sx q[0];
rz(-1.5710314) q[0];
sx q[0];
rz(-1.4192441) q[0];
x q[1];
rz(-1.542369) q[2];
sx q[2];
rz(-2.5027983) q[2];
sx q[2];
rz(-0.5109238) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6533493) q[1];
sx q[1];
rz(-2.0470139) q[1];
sx q[1];
rz(-0.32250065) q[1];
rz(-1.7436036) q[3];
sx q[3];
rz(-0.51522845) q[3];
sx q[3];
rz(-2.990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8511054) q[2];
sx q[2];
rz(-1.9600211) q[2];
sx q[2];
rz(0.30430749) q[2];
rz(2.815222) q[3];
sx q[3];
rz(-0.54309741) q[3];
sx q[3];
rz(-0.91199005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07311634) q[0];
sx q[0];
rz(-2.731972) q[0];
sx q[0];
rz(2.2090744) q[0];
rz(0.069123507) q[1];
sx q[1];
rz(-2.9239475) q[1];
sx q[1];
rz(2.1014012) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.067834082) q[0];
sx q[0];
rz(-2.375575) q[0];
sx q[0];
rz(0.037184663) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1625453) q[2];
sx q[2];
rz(-0.30610105) q[2];
sx q[2];
rz(3.0476168) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5835147) q[1];
sx q[1];
rz(-0.78418523) q[1];
sx q[1];
rz(0.23434831) q[1];
rz(-pi) q[2];
rz(0.64994855) q[3];
sx q[3];
rz(-0.81226617) q[3];
sx q[3];
rz(2.5358729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0065464) q[2];
sx q[2];
rz(-2.4032205) q[2];
sx q[2];
rz(-0.73299232) q[2];
rz(1.0829571) q[3];
sx q[3];
rz(-1.0561918) q[3];
sx q[3];
rz(-1.0068896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5873544) q[0];
sx q[0];
rz(-0.08389689) q[0];
sx q[0];
rz(0.49302897) q[0];
rz(-2.9250277) q[1];
sx q[1];
rz(-1.1410057) q[1];
sx q[1];
rz(2.1176178) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42070779) q[0];
sx q[0];
rz(-2.5455089) q[0];
sx q[0];
rz(-3.1203624) q[0];
x q[1];
rz(-0.31957133) q[2];
sx q[2];
rz(-0.34863483) q[2];
sx q[2];
rz(0.92952585) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1223381) q[1];
sx q[1];
rz(-2.755909) q[1];
sx q[1];
rz(-2.7606439) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8423457) q[3];
sx q[3];
rz(-0.95029921) q[3];
sx q[3];
rz(-2.3503691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0563125) q[2];
sx q[2];
rz(-1.5093426) q[2];
sx q[2];
rz(2.4658266) q[2];
rz(2.7285649) q[3];
sx q[3];
rz(-1.4512117) q[3];
sx q[3];
rz(-0.54615027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6290879) q[0];
sx q[0];
rz(-0.65852037) q[0];
sx q[0];
rz(0.034544695) q[0];
rz(-0.054556219) q[1];
sx q[1];
rz(-1.9787534) q[1];
sx q[1];
rz(0.82130718) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1217348) q[0];
sx q[0];
rz(-0.99911753) q[0];
sx q[0];
rz(-1.4522533) q[0];
rz(-pi) q[1];
rz(-0.088772687) q[2];
sx q[2];
rz(-1.5615511) q[2];
sx q[2];
rz(-2.451596) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8470314) q[1];
sx q[1];
rz(-2.9518513) q[1];
sx q[1];
rz(2.4824597) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6590632) q[3];
sx q[3];
rz(-2.165613) q[3];
sx q[3];
rz(1.1728668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36046946) q[2];
sx q[2];
rz(-1.0162153) q[2];
sx q[2];
rz(-0.13295573) q[2];
rz(-0.41108701) q[3];
sx q[3];
rz(-1.5186331) q[3];
sx q[3];
rz(-2.0558004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046831176) q[0];
sx q[0];
rz(-0.60698858) q[0];
sx q[0];
rz(1.2588311) q[0];
rz(-1.6350485) q[1];
sx q[1];
rz(-2.4848487) q[1];
sx q[1];
rz(-2.3283995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98018089) q[0];
sx q[0];
rz(-1.6488355) q[0];
sx q[0];
rz(1.6059931) q[0];
x q[1];
rz(1.7178463) q[2];
sx q[2];
rz(-2.4644901) q[2];
sx q[2];
rz(-0.33237095) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30505667) q[1];
sx q[1];
rz(-2.146135) q[1];
sx q[1];
rz(2.1161377) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4370632) q[3];
sx q[3];
rz(-2.0195035) q[3];
sx q[3];
rz(0.98380849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4930111) q[2];
sx q[2];
rz(-2.0040671) q[2];
sx q[2];
rz(-2.0900334) q[2];
rz(1.018853) q[3];
sx q[3];
rz(-2.2994883) q[3];
sx q[3];
rz(-0.65783182) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9857585) q[0];
sx q[0];
rz(-0.65407615) q[0];
sx q[0];
rz(0.88055897) q[0];
rz(-1.8219933) q[1];
sx q[1];
rz(-1.1450014) q[1];
sx q[1];
rz(-3.0897279) q[1];
rz(2.8116799) q[2];
sx q[2];
rz(-2.920426) q[2];
sx q[2];
rz(-2.5830808) q[2];
rz(1.5108091) q[3];
sx q[3];
rz(-1.7675478) q[3];
sx q[3];
rz(0.02789733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
