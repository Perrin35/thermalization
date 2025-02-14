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
rz(-1.2404233) q[0];
sx q[0];
rz(-1.197553) q[0];
sx q[0];
rz(-0.21633202) q[0];
rz(4.0989838) q[1];
sx q[1];
rz(5.6070072) q[1];
sx q[1];
rz(11.054872) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51115655) q[0];
sx q[0];
rz(-1.5611708) q[0];
sx q[0];
rz(-2.5853392) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6681603) q[2];
sx q[2];
rz(-1.8634999) q[2];
sx q[2];
rz(-2.6930489) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.87147698) q[1];
sx q[1];
rz(-1.2387215) q[1];
sx q[1];
rz(-2.6144652) q[1];
rz(0.39325289) q[3];
sx q[3];
rz(-2.234708) q[3];
sx q[3];
rz(0.25379405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.034885255) q[2];
sx q[2];
rz(-0.35197508) q[2];
sx q[2];
rz(1.0484288) q[2];
rz(2.9599221) q[3];
sx q[3];
rz(-2.1766267) q[3];
sx q[3];
rz(0.65339965) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6492017) q[0];
sx q[0];
rz(-0.94434706) q[0];
sx q[0];
rz(0.44678584) q[0];
rz(-0.88042879) q[1];
sx q[1];
rz(-1.7767521) q[1];
sx q[1];
rz(0.78278881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.114417) q[0];
sx q[0];
rz(-1.3265298) q[0];
sx q[0];
rz(3.0784803) q[0];
rz(-pi) q[1];
rz(2.5249285) q[2];
sx q[2];
rz(-1.2149802) q[2];
sx q[2];
rz(2.7105376) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.782541) q[1];
sx q[1];
rz(-0.56582574) q[1];
sx q[1];
rz(0.58360175) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47329013) q[3];
sx q[3];
rz(-2.1968699) q[3];
sx q[3];
rz(-1.9260709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.11144) q[2];
sx q[2];
rz(-1.4806662) q[2];
sx q[2];
rz(-2.1622369) q[2];
rz(2.7566946) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(2.9773007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48591831) q[0];
sx q[0];
rz(-0.05412183) q[0];
sx q[0];
rz(0.78980494) q[0];
rz(-2.9549331) q[1];
sx q[1];
rz(-1.7402382) q[1];
sx q[1];
rz(-0.99367118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2659205) q[0];
sx q[0];
rz(-1.3215995) q[0];
sx q[0];
rz(-2.5795291) q[0];
x q[1];
rz(-2.5285401) q[2];
sx q[2];
rz(-1.7568577) q[2];
sx q[2];
rz(1.8579872) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4504272) q[1];
sx q[1];
rz(-1.9948927) q[1];
sx q[1];
rz(2.0753839) q[1];
x q[2];
rz(0.72088269) q[3];
sx q[3];
rz(-1.7768806) q[3];
sx q[3];
rz(-0.64180798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.32187244) q[2];
sx q[2];
rz(-2.0237782) q[2];
sx q[2];
rz(2.7703088) q[2];
rz(-2.7927981) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(-2.546052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45469859) q[0];
sx q[0];
rz(-2.1214387) q[0];
sx q[0];
rz(1.4042847) q[0];
rz(2.4644201) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(-1.4926532) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9335404) q[0];
sx q[0];
rz(-1.786199) q[0];
sx q[0];
rz(0.24788863) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8437496) q[2];
sx q[2];
rz(-1.1561511) q[2];
sx q[2];
rz(2.9306102) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4863534) q[1];
sx q[1];
rz(-0.35105536) q[1];
sx q[1];
rz(0.69610657) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31217137) q[3];
sx q[3];
rz(-2.2709284) q[3];
sx q[3];
rz(-0.90538245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6877785) q[2];
sx q[2];
rz(-0.9684338) q[2];
sx q[2];
rz(-2.7933534) q[2];
rz(-1.4549152) q[3];
sx q[3];
rz(-1.7125407) q[3];
sx q[3];
rz(1.3454364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1763879) q[0];
sx q[0];
rz(-2.5768953) q[0];
sx q[0];
rz(-2.2494466) q[0];
rz(0.46420321) q[1];
sx q[1];
rz(-1.2449539) q[1];
sx q[1];
rz(-1.2947327) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1682353) q[0];
sx q[0];
rz(-1.1001612) q[0];
sx q[0];
rz(-0.55460978) q[0];
x q[1];
rz(-2.928435) q[2];
sx q[2];
rz(-2.5937754) q[2];
sx q[2];
rz(2.167706) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0820992) q[1];
sx q[1];
rz(-1.1319185) q[1];
sx q[1];
rz(2.782269) q[1];
x q[2];
rz(-0.89035676) q[3];
sx q[3];
rz(-2.0687752) q[3];
sx q[3];
rz(0.80163664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3225473) q[2];
sx q[2];
rz(-2.5307541) q[2];
sx q[2];
rz(-0.56582212) q[2];
rz(-3.0602509) q[3];
sx q[3];
rz(-2.1606725) q[3];
sx q[3];
rz(2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.464798) q[0];
sx q[0];
rz(-3.0749574) q[0];
sx q[0];
rz(-1.5860522) q[0];
rz(2.0809035) q[1];
sx q[1];
rz(-1.5763177) q[1];
sx q[1];
rz(-0.63180822) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.466296) q[0];
sx q[0];
rz(-1.491437) q[0];
sx q[0];
rz(2.9064889) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6423694) q[2];
sx q[2];
rz(-0.92485917) q[2];
sx q[2];
rz(-1.8725841) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8579156) q[1];
sx q[1];
rz(-1.2322958) q[1];
sx q[1];
rz(1.0033016) q[1];
rz(2.4580509) q[3];
sx q[3];
rz(-1.6850796) q[3];
sx q[3];
rz(-1.5147234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1598728) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(-2.6427606) q[2];
rz(-1.8286797) q[3];
sx q[3];
rz(-0.73176089) q[3];
sx q[3];
rz(-1.7074728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.53443921) q[0];
sx q[0];
rz(-2.7767015) q[0];
sx q[0];
rz(-2.4420807) q[0];
rz(2.6761159) q[1];
sx q[1];
rz(-0.87124467) q[1];
sx q[1];
rz(0.93719283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95109601) q[0];
sx q[0];
rz(-1.3590004) q[0];
sx q[0];
rz(-1.3986716) q[0];
rz(-pi) q[1];
rz(2.7363051) q[2];
sx q[2];
rz(-2.8817085) q[2];
sx q[2];
rz(2.0854307) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4027462) q[1];
sx q[1];
rz(-0.35772309) q[1];
sx q[1];
rz(2.2132232) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8392467) q[3];
sx q[3];
rz(-1.8519326) q[3];
sx q[3];
rz(-2.9535563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.84404868) q[2];
sx q[2];
rz(-1.301845) q[2];
sx q[2];
rz(0.37332264) q[2];
rz(1.1897872) q[3];
sx q[3];
rz(-0.50896421) q[3];
sx q[3];
rz(0.51026195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74476403) q[0];
sx q[0];
rz(-2.0592392) q[0];
sx q[0];
rz(2.3936791) q[0];
rz(0.76639908) q[1];
sx q[1];
rz(-2.8728569) q[1];
sx q[1];
rz(-3.1386197) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.317694) q[0];
sx q[0];
rz(-0.71461535) q[0];
sx q[0];
rz(-1.8665642) q[0];
rz(1.4998798) q[2];
sx q[2];
rz(-1.5902963) q[2];
sx q[2];
rz(1.2649346) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6312478) q[1];
sx q[1];
rz(-1.6142134) q[1];
sx q[1];
rz(1.4166635) q[1];
rz(-pi) q[2];
x q[2];
rz(0.19488867) q[3];
sx q[3];
rz(-1.7656544) q[3];
sx q[3];
rz(-3.0611567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.030674) q[2];
sx q[2];
rz(-1.3838394) q[2];
sx q[2];
rz(-1.0670916) q[2];
rz(3.057632) q[3];
sx q[3];
rz(-0.48560086) q[3];
sx q[3];
rz(0.77897227) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79754168) q[0];
sx q[0];
rz(-2.2028956) q[0];
sx q[0];
rz(0.10636605) q[0];
rz(2.1616409) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(-0.76593691) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2666616) q[0];
sx q[0];
rz(-0.70737544) q[0];
sx q[0];
rz(-2.4445663) q[0];
x q[1];
rz(-2.9418403) q[2];
sx q[2];
rz(-1.6589266) q[2];
sx q[2];
rz(-0.55921184) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76607219) q[1];
sx q[1];
rz(-0.67691708) q[1];
sx q[1];
rz(3.0860597) q[1];
rz(-0.25772734) q[3];
sx q[3];
rz(-0.014583909) q[3];
sx q[3];
rz(-2.7590883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6672259) q[2];
sx q[2];
rz(-2.6211278) q[2];
sx q[2];
rz(-0.70029798) q[2];
rz(0.66655603) q[3];
sx q[3];
rz(-1.3349814) q[3];
sx q[3];
rz(2.0535645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661082) q[0];
sx q[0];
rz(-2.1145144) q[0];
sx q[0];
rz(-1.3960557) q[0];
rz(1.4736157) q[1];
sx q[1];
rz(-1.2584078) q[1];
sx q[1];
rz(2.4748763) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7233492) q[0];
sx q[0];
rz(-2.1921792) q[0];
sx q[0];
rz(-2.3410812) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.780203) q[2];
sx q[2];
rz(-0.69212428) q[2];
sx q[2];
rz(-0.036594242) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8620876) q[1];
sx q[1];
rz(-2.7976329) q[1];
sx q[1];
rz(0.39744795) q[1];
rz(0.57228831) q[3];
sx q[3];
rz(-1.0977355) q[3];
sx q[3];
rz(0.285404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.066976808) q[2];
sx q[2];
rz(-0.65698996) q[2];
sx q[2];
rz(-2.2508049) q[2];
rz(-2.7095419) q[3];
sx q[3];
rz(-2.0712349) q[3];
sx q[3];
rz(0.74705684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4074832) q[0];
sx q[0];
rz(-1.4984087) q[0];
sx q[0];
rz(1.8468504) q[0];
rz(-1.8941849) q[1];
sx q[1];
rz(-0.71949646) q[1];
sx q[1];
rz(1.234642) q[1];
rz(0.16130372) q[2];
sx q[2];
rz(-1.9111173) q[2];
sx q[2];
rz(2.5572122) q[2];
rz(2.6525146) q[3];
sx q[3];
rz(-1.6001971) q[3];
sx q[3];
rz(-1.8501545) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
