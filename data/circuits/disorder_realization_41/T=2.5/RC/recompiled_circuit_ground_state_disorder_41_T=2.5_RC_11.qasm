OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.95076686) q[0];
sx q[0];
rz(-1.8523676) q[0];
sx q[0];
rz(-3.0486795) q[0];
rz(-0.095280401) q[1];
sx q[1];
rz(-0.73268259) q[1];
sx q[1];
rz(1.9021775) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026442083) q[0];
sx q[0];
rz(-1.7998994) q[0];
sx q[0];
rz(1.312448) q[0];
x q[1];
rz(-2.8777166) q[2];
sx q[2];
rz(-1.3812764) q[2];
sx q[2];
rz(-2.5352728) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3500786) q[1];
sx q[1];
rz(-2.8100612) q[1];
sx q[1];
rz(3.0505807) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4437583) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(-2.7891911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4702845) q[2];
sx q[2];
rz(-1.6613864) q[2];
sx q[2];
rz(1.3489464) q[2];
rz(0.74364439) q[3];
sx q[3];
rz(-0.27739224) q[3];
sx q[3];
rz(-0.17717895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.141356) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(0.29362383) q[0];
rz(2.1108421) q[1];
sx q[1];
rz(-1.3615969) q[1];
sx q[1];
rz(2.7986599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2156649) q[0];
sx q[0];
rz(-0.56680381) q[0];
sx q[0];
rz(1.180992) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9428217) q[2];
sx q[2];
rz(-1.6918394) q[2];
sx q[2];
rz(-2.4189359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56655069) q[1];
sx q[1];
rz(-2.2600265) q[1];
sx q[1];
rz(1.291996) q[1];
rz(-pi) q[2];
rz(2.2003056) q[3];
sx q[3];
rz(-1.7751179) q[3];
sx q[3];
rz(2.2245537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.12884101) q[2];
sx q[2];
rz(-1.3920471) q[2];
sx q[2];
rz(2.3823605) q[2];
rz(-0.020708474) q[3];
sx q[3];
rz(-1.8755251) q[3];
sx q[3];
rz(2.7032963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6193806) q[0];
sx q[0];
rz(-0.601957) q[0];
sx q[0];
rz(0.0044599175) q[0];
rz(2.4941173) q[1];
sx q[1];
rz(-0.59440333) q[1];
sx q[1];
rz(-0.78537816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42703907) q[0];
sx q[0];
rz(-1.4903965) q[0];
sx q[0];
rz(-2.9820082) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0702281) q[2];
sx q[2];
rz(-2.7874261) q[2];
sx q[2];
rz(-0.11812299) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.33357802) q[1];
sx q[1];
rz(-1.6393264) q[1];
sx q[1];
rz(-1.0295111) q[1];
rz(2.3397869) q[3];
sx q[3];
rz(-0.84322819) q[3];
sx q[3];
rz(1.8319195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.37041) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(-1.3578337) q[2];
rz(-2.8010098) q[3];
sx q[3];
rz(-0.95652306) q[3];
sx q[3];
rz(-2.3497439) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1444645) q[0];
sx q[0];
rz(-1.0839394) q[0];
sx q[0];
rz(2.2564364) q[0];
rz(-0.74686933) q[1];
sx q[1];
rz(-2.5179458) q[1];
sx q[1];
rz(1.0964099) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089450739) q[0];
sx q[0];
rz(-0.47596395) q[0];
sx q[0];
rz(-2.3030998) q[0];
rz(-pi) q[1];
rz(-1.5113513) q[2];
sx q[2];
rz(-1.6013833) q[2];
sx q[2];
rz(2.2074043) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.29624048) q[1];
sx q[1];
rz(-1.3120323) q[1];
sx q[1];
rz(-2.6654763) q[1];
rz(-pi) q[2];
rz(-1.0616333) q[3];
sx q[3];
rz(-0.90583767) q[3];
sx q[3];
rz(0.10275118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61802822) q[2];
sx q[2];
rz(-2.9310493) q[2];
sx q[2];
rz(-0.024624126) q[2];
rz(0.63557449) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(2.2583101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4516975) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(0.46689335) q[0];
rz(3.0310071) q[1];
sx q[1];
rz(-2.6778335) q[1];
sx q[1];
rz(-1.0708403) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3239646) q[0];
sx q[0];
rz(-1.1187129) q[0];
sx q[0];
rz(-0.26497193) q[0];
rz(-1.3883287) q[2];
sx q[2];
rz(-1.8844205) q[2];
sx q[2];
rz(-0.30280606) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3339798) q[1];
sx q[1];
rz(-2.4845504) q[1];
sx q[1];
rz(0.64589898) q[1];
rz(0.73770817) q[3];
sx q[3];
rz(-2.2900351) q[3];
sx q[3];
rz(-1.4054474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0232627) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(0.85232097) q[2];
rz(0.037633745) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(2.5467303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1325876) q[0];
sx q[0];
rz(-1.9629033) q[0];
sx q[0];
rz(1.2888541) q[0];
rz(1.5124849) q[1];
sx q[1];
rz(-1.1237203) q[1];
sx q[1];
rz(-1.1423133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4415671) q[0];
sx q[0];
rz(-1.4730994) q[0];
sx q[0];
rz(2.7872681) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3630168) q[2];
sx q[2];
rz(-1.6315952) q[2];
sx q[2];
rz(-0.5375934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2494275) q[1];
sx q[1];
rz(-1.2293503) q[1];
sx q[1];
rz(-1.83431) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31123881) q[3];
sx q[3];
rz(-1.4025619) q[3];
sx q[3];
rz(2.5995035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3219354) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(3.0211871) q[2];
rz(1.1558007) q[3];
sx q[3];
rz(-1.5294231) q[3];
sx q[3];
rz(0.22404484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33893809) q[0];
sx q[0];
rz(-1.3405565) q[0];
sx q[0];
rz(-1.6876203) q[0];
rz(-1.5974207) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(-2.6905751) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16674834) q[0];
sx q[0];
rz(-1.5112425) q[0];
sx q[0];
rz(-0.3015428) q[0];
rz(-pi) q[1];
rz(1.0181997) q[2];
sx q[2];
rz(-2.7180053) q[2];
sx q[2];
rz(0.46679631) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1731733) q[1];
sx q[1];
rz(-2.0035158) q[1];
sx q[1];
rz(-1.7041901) q[1];
rz(-pi) q[2];
rz(-3.1322828) q[3];
sx q[3];
rz(-2.3311989) q[3];
sx q[3];
rz(-1.1610462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94577998) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(-1.7808524) q[2];
rz(2.1360548) q[3];
sx q[3];
rz(-1.9660299) q[3];
sx q[3];
rz(1.3895234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239531) q[0];
sx q[0];
rz(-3.0160976) q[0];
sx q[0];
rz(0.70607591) q[0];
rz(-1.5438682) q[1];
sx q[1];
rz(-1.7192625) q[1];
sx q[1];
rz(-2.2854038) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5917011) q[0];
sx q[0];
rz(-0.020792637) q[0];
sx q[0];
rz(-1.2077622) q[0];
rz(-pi) q[1];
rz(-2.0308308) q[2];
sx q[2];
rz(-0.80321124) q[2];
sx q[2];
rz(0.75424657) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4256712) q[1];
sx q[1];
rz(-1.7737264) q[1];
sx q[1];
rz(-2.5589866) q[1];
x q[2];
rz(2.059465) q[3];
sx q[3];
rz(-1.6889945) q[3];
sx q[3];
rz(1.4335872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.84264821) q[2];
sx q[2];
rz(-1.7405258) q[2];
sx q[2];
rz(-2.974158) q[2];
rz(1.5873448) q[3];
sx q[3];
rz(-2.3685679) q[3];
sx q[3];
rz(2.1412444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7635968) q[0];
sx q[0];
rz(-1.6908228) q[0];
sx q[0];
rz(2.7698621) q[0];
rz(1.7077712) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(-2.8415714) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094291501) q[0];
sx q[0];
rz(-0.30891788) q[0];
sx q[0];
rz(-2.7500217) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5302363) q[2];
sx q[2];
rz(-2.1475361) q[2];
sx q[2];
rz(-0.64740136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3206818) q[1];
sx q[1];
rz(-1.4467941) q[1];
sx q[1];
rz(2.4895666) q[1];
rz(-pi) q[2];
rz(-2.2737502) q[3];
sx q[3];
rz(-1.4757089) q[3];
sx q[3];
rz(-0.33012182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.47664777) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(2.7139968) q[2];
rz(-1.556373) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(-0.75470406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72117358) q[0];
sx q[0];
rz(-2.6009646) q[0];
sx q[0];
rz(2.6721201) q[0];
rz(2.2162614) q[1];
sx q[1];
rz(-1.0877437) q[1];
sx q[1];
rz(3.1184149) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4172702) q[0];
sx q[0];
rz(-2.5886154) q[0];
sx q[0];
rz(0.31674762) q[0];
rz(-pi) q[1];
rz(-2.2397348) q[2];
sx q[2];
rz(-2.0624071) q[2];
sx q[2];
rz(-0.41504809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7851023) q[1];
sx q[1];
rz(-1.7447326) q[1];
sx q[1];
rz(-2.9746303) q[1];
x q[2];
rz(-1.7816824) q[3];
sx q[3];
rz(-2.4004705) q[3];
sx q[3];
rz(-1.1531342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2716486) q[2];
sx q[2];
rz(-1.9359438) q[2];
sx q[2];
rz(2.6507157) q[2];
rz(1.0902181) q[3];
sx q[3];
rz(-0.96929437) q[3];
sx q[3];
rz(0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8563817) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(2.8648227) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(-0.3666975) q[2];
sx q[2];
rz(-0.72281217) q[2];
sx q[2];
rz(-2.6642961) q[2];
rz(-0.53448813) q[3];
sx q[3];
rz(-1.5535003) q[3];
sx q[3];
rz(1.5737892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
