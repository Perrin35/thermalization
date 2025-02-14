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
rz(-0.99875206) q[0];
sx q[0];
rz(-0.771703) q[0];
sx q[0];
rz(1.9089215) q[0];
rz(-1.9219037) q[1];
sx q[1];
rz(-0.86644679) q[1];
sx q[1];
rz(-2.7629857) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50921665) q[0];
sx q[0];
rz(-1.6086744) q[0];
sx q[0];
rz(1.3404472) q[0];
x q[1];
rz(-0.98388715) q[2];
sx q[2];
rz(-0.68418938) q[2];
sx q[2];
rz(0.049287576) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9512082) q[1];
sx q[1];
rz(-2.4255014) q[1];
sx q[1];
rz(-0.33204099) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55128184) q[3];
sx q[3];
rz(-2.4935185) q[3];
sx q[3];
rz(0.66873811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1401225) q[2];
sx q[2];
rz(-1.9688316) q[2];
sx q[2];
rz(2.6918461) q[2];
rz(-2.5887515) q[3];
sx q[3];
rz(-0.55838412) q[3];
sx q[3];
rz(0.21137485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1028035) q[0];
sx q[0];
rz(-1.1454104) q[0];
sx q[0];
rz(0.75622028) q[0];
rz(-2.4839632) q[1];
sx q[1];
rz(-0.88784528) q[1];
sx q[1];
rz(2.5530946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5969944) q[0];
sx q[0];
rz(-2.2203494) q[0];
sx q[0];
rz(-2.0429595) q[0];
rz(-pi) q[1];
rz(-1.8075077) q[2];
sx q[2];
rz(-0.76067096) q[2];
sx q[2];
rz(-2.0784645) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64526833) q[1];
sx q[1];
rz(-0.2926188) q[1];
sx q[1];
rz(-1.5436548) q[1];
rz(-pi) q[2];
rz(1.9733834) q[3];
sx q[3];
rz(-0.63791554) q[3];
sx q[3];
rz(2.0762553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6077891) q[2];
sx q[2];
rz(-1.6705931) q[2];
sx q[2];
rz(3.0652453) q[2];
rz(-0.40147436) q[3];
sx q[3];
rz(-2.6186826) q[3];
sx q[3];
rz(2.0126066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56979316) q[0];
sx q[0];
rz(-2.5558668) q[0];
sx q[0];
rz(0.3366003) q[0];
rz(-2.1062689) q[1];
sx q[1];
rz(-0.88408771) q[1];
sx q[1];
rz(-3.091277) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53165689) q[0];
sx q[0];
rz(-0.7748403) q[0];
sx q[0];
rz(0.70506238) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5467954) q[2];
sx q[2];
rz(-1.7414879) q[2];
sx q[2];
rz(-1.6015944) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2256098) q[1];
sx q[1];
rz(-0.40357737) q[1];
sx q[1];
rz(-0.55761375) q[1];
x q[2];
rz(0.21304275) q[3];
sx q[3];
rz(-2.1318949) q[3];
sx q[3];
rz(1.9764501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7670373) q[2];
sx q[2];
rz(-0.87696806) q[2];
sx q[2];
rz(0.98131895) q[2];
rz(2.0902925) q[3];
sx q[3];
rz(-1.4562573) q[3];
sx q[3];
rz(0.013896996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0392847) q[0];
sx q[0];
rz(-1.4816875) q[0];
sx q[0];
rz(-2.1639977) q[0];
rz(-0.54999894) q[1];
sx q[1];
rz(-2.1606052) q[1];
sx q[1];
rz(-0.94211334) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4956168) q[0];
sx q[0];
rz(-1.2496523) q[0];
sx q[0];
rz(0.34713388) q[0];
rz(-0.35772985) q[2];
sx q[2];
rz(-1.5204645) q[2];
sx q[2];
rz(0.77046399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4244183) q[1];
sx q[1];
rz(-2.5150329) q[1];
sx q[1];
rz(-1.8612629) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1051543) q[3];
sx q[3];
rz(-2.0490408) q[3];
sx q[3];
rz(-2.382232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1407239) q[2];
sx q[2];
rz(-2.7694747) q[2];
sx q[2];
rz(-1.6288527) q[2];
rz(0.77423972) q[3];
sx q[3];
rz(-2.1858678) q[3];
sx q[3];
rz(0.16404185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0999488) q[0];
sx q[0];
rz(-2.994717) q[0];
sx q[0];
rz(-0.80648333) q[0];
rz(2.3355314) q[1];
sx q[1];
rz(-1.1763923) q[1];
sx q[1];
rz(2.3890266) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0092750015) q[0];
sx q[0];
rz(-0.89269367) q[0];
sx q[0];
rz(0.14831774) q[0];
rz(-1.5335778) q[2];
sx q[2];
rz(-1.9646837) q[2];
sx q[2];
rz(-1.5348127) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6347152) q[1];
sx q[1];
rz(-1.1487097) q[1];
sx q[1];
rz(-1.6812801) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36461501) q[3];
sx q[3];
rz(-2.7066253) q[3];
sx q[3];
rz(1.9221969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8807184) q[2];
sx q[2];
rz(-1.4148834) q[2];
sx q[2];
rz(0.073337642) q[2];
rz(0.65555769) q[3];
sx q[3];
rz(-0.95584241) q[3];
sx q[3];
rz(-0.59371322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23621479) q[0];
sx q[0];
rz(-2.0423934) q[0];
sx q[0];
rz(2.4429876) q[0];
rz(-1.3911432) q[1];
sx q[1];
rz(-0.51934424) q[1];
sx q[1];
rz(-3.0379675) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9148283) q[0];
sx q[0];
rz(-1.6092759) q[0];
sx q[0];
rz(1.5917042) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.03662911) q[2];
sx q[2];
rz(-2.0994791) q[2];
sx q[2];
rz(1.5122459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0558529) q[1];
sx q[1];
rz(-0.21634783) q[1];
sx q[1];
rz(2.0592709) q[1];
rz(1.7781939) q[3];
sx q[3];
rz(-1.816664) q[3];
sx q[3];
rz(1.8519925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.708272) q[2];
sx q[2];
rz(-2.2430113) q[2];
sx q[2];
rz(1.9898604) q[2];
rz(-0.084130675) q[3];
sx q[3];
rz(-0.79456544) q[3];
sx q[3];
rz(2.7860723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5750835) q[0];
sx q[0];
rz(-2.7783448) q[0];
sx q[0];
rz(1.890924) q[0];
rz(-1.7364712) q[1];
sx q[1];
rz(-0.57650081) q[1];
sx q[1];
rz(2.0723431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.273928) q[0];
sx q[0];
rz(-1.3153512) q[0];
sx q[0];
rz(2.5219265) q[0];
rz(-pi) q[1];
rz(-1.5450412) q[2];
sx q[2];
rz(-0.47821486) q[2];
sx q[2];
rz(1.1104294) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7916966) q[1];
sx q[1];
rz(-2.2772067) q[1];
sx q[1];
rz(3.033925) q[1];
x q[2];
rz(0.077344374) q[3];
sx q[3];
rz(-1.8048022) q[3];
sx q[3];
rz(1.6868433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8413267) q[2];
sx q[2];
rz(-1.0043251) q[2];
sx q[2];
rz(-2.0036073) q[2];
rz(0.56002069) q[3];
sx q[3];
rz(-1.0617826) q[3];
sx q[3];
rz(-1.2573857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1244125) q[0];
sx q[0];
rz(-0.37473285) q[0];
sx q[0];
rz(2.4770233) q[0];
rz(-2.7542704) q[1];
sx q[1];
rz(-2.1699984) q[1];
sx q[1];
rz(1.1727146) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5636503) q[0];
sx q[0];
rz(-1.1563753) q[0];
sx q[0];
rz(2.2060288) q[0];
rz(-pi) q[1];
rz(-2.6792999) q[2];
sx q[2];
rz(-1.3707118) q[2];
sx q[2];
rz(-2.878396) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.203158) q[1];
sx q[1];
rz(-2.3745792) q[1];
sx q[1];
rz(-2.2545283) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6608638) q[3];
sx q[3];
rz(-0.58739118) q[3];
sx q[3];
rz(1.4241854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6465801) q[2];
sx q[2];
rz(-2.6436372) q[2];
sx q[2];
rz(-0.23769561) q[2];
rz(-0.85463917) q[3];
sx q[3];
rz(-1.6228638) q[3];
sx q[3];
rz(2.2304227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89331996) q[0];
sx q[0];
rz(-1.3398291) q[0];
sx q[0];
rz(2.7082537) q[0];
rz(1.8386819) q[1];
sx q[1];
rz(-1.358526) q[1];
sx q[1];
rz(-0.059159577) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3179037) q[0];
sx q[0];
rz(-1.5971703) q[0];
sx q[0];
rz(-1.8120873) q[0];
rz(-pi) q[1];
rz(-1.3174345) q[2];
sx q[2];
rz(-2.4640016) q[2];
sx q[2];
rz(2.3357725) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6530919) q[1];
sx q[1];
rz(-1.8951804) q[1];
sx q[1];
rz(-1.2725194) q[1];
rz(1.9800277) q[3];
sx q[3];
rz(-0.61823119) q[3];
sx q[3];
rz(-0.78200862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9383508) q[2];
sx q[2];
rz(-1.4138736) q[2];
sx q[2];
rz(2.4436277) q[2];
rz(-0.62567726) q[3];
sx q[3];
rz(-2.8901849) q[3];
sx q[3];
rz(-1.8184557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0980001) q[0];
sx q[0];
rz(-2.5830144) q[0];
sx q[0];
rz(0.25398764) q[0];
rz(0.99098539) q[1];
sx q[1];
rz(-1.5373983) q[1];
sx q[1];
rz(3.1414247) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5714665) q[0];
sx q[0];
rz(-1.1227221) q[0];
sx q[0];
rz(-2.6506533) q[0];
x q[1];
rz(-1.1824047) q[2];
sx q[2];
rz(-0.12586181) q[2];
sx q[2];
rz(-2.2110155) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.907577) q[1];
sx q[1];
rz(-1.8157417) q[1];
sx q[1];
rz(2.8744389) q[1];
rz(0.17642085) q[3];
sx q[3];
rz(-1.9690445) q[3];
sx q[3];
rz(-0.12052025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1508472) q[2];
sx q[2];
rz(-1.3517697) q[2];
sx q[2];
rz(-0.051648971) q[2];
rz(2.1150186) q[3];
sx q[3];
rz(-2.5995422) q[3];
sx q[3];
rz(-0.76977229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50292618) q[0];
sx q[0];
rz(-2.2421056) q[0];
sx q[0];
rz(-2.5594287) q[0];
rz(0.40099405) q[1];
sx q[1];
rz(-2.4759226) q[1];
sx q[1];
rz(-0.84025875) q[1];
rz(2.6225093) q[2];
sx q[2];
rz(-1.097403) q[2];
sx q[2];
rz(2.335142) q[2];
rz(-2.9874728) q[3];
sx q[3];
rz(-1.6748292) q[3];
sx q[3];
rz(-2.0451121) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
