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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5372972) q[0];
sx q[0];
rz(-1.3193466) q[0];
sx q[0];
rz(-0.23668134) q[0];
rz(-1.3746512) q[2];
sx q[2];
rz(-1.3117547) q[2];
sx q[2];
rz(-2.1262622) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.8653633) q[1];
sx q[1];
rz(-1.5412093) q[1];
sx q[1];
rz(-2.8113356) q[1];
x q[2];
rz(2.4437583) q[3];
sx q[3];
rz(-2.4437332) q[3];
sx q[3];
rz(-2.7891911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4702845) q[2];
sx q[2];
rz(-1.4802063) q[2];
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
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0002366) q[0];
sx q[0];
rz(-1.5970705) q[0];
sx q[0];
rz(-0.29362383) q[0];
rz(1.0307505) q[1];
sx q[1];
rz(-1.3615969) q[1];
sx q[1];
rz(-2.7986599) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4727508) q[0];
sx q[0];
rz(-2.0905417) q[0];
sx q[0];
rz(2.9042871) q[0];
rz(-0.12983506) q[2];
sx q[2];
rz(-1.2016244) q[2];
sx q[2];
rz(-0.8952199) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.82417497) q[1];
sx q[1];
rz(-1.7848099) q[1];
sx q[1];
rz(0.70862464) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2003056) q[3];
sx q[3];
rz(-1.7751179) q[3];
sx q[3];
rz(-2.2245537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12884101) q[2];
sx q[2];
rz(-1.7495456) q[2];
sx q[2];
rz(0.7592321) q[2];
rz(0.020708474) q[3];
sx q[3];
rz(-1.2660675) q[3];
sx q[3];
rz(2.7032963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.52221209) q[0];
sx q[0];
rz(-0.601957) q[0];
sx q[0];
rz(0.0044599175) q[0];
rz(-2.4941173) q[1];
sx q[1];
rz(-2.5471893) q[1];
sx q[1];
rz(2.3562145) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1566832) q[0];
sx q[0];
rz(-1.4117318) q[0];
sx q[0];
rz(-1.6522264) q[0];
rz(-pi) q[1];
rz(1.8846604) q[2];
sx q[2];
rz(-1.7376657) q[2];
sx q[2];
rz(2.1618869) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.33357802) q[1];
sx q[1];
rz(-1.6393264) q[1];
sx q[1];
rz(2.1120815) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80180577) q[3];
sx q[3];
rz(-2.2983645) q[3];
sx q[3];
rz(1.8319195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.37041) q[2];
sx q[2];
rz(-2.2213171) q[2];
sx q[2];
rz(1.7837589) q[2];
rz(-2.8010098) q[3];
sx q[3];
rz(-2.1850696) q[3];
sx q[3];
rz(-0.79184872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1444645) q[0];
sx q[0];
rz(-2.0576532) q[0];
sx q[0];
rz(0.88515627) q[0];
rz(-2.3947233) q[1];
sx q[1];
rz(-0.62364686) q[1];
sx q[1];
rz(1.0964099) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0521419) q[0];
sx q[0];
rz(-2.6656287) q[0];
sx q[0];
rz(0.83849283) q[0];
x q[1];
rz(1.5113513) q[2];
sx q[2];
rz(-1.6013833) q[2];
sx q[2];
rz(-2.2074043) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1433574) q[1];
sx q[1];
rz(-1.1117742) q[1];
sx q[1];
rz(1.281339) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55620749) q[3];
sx q[3];
rz(-0.81327754) q[3];
sx q[3];
rz(2.5089056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5235644) q[2];
sx q[2];
rz(-0.21054331) q[2];
sx q[2];
rz(-3.1169685) q[2];
rz(0.63557449) q[3];
sx q[3];
rz(-2.1533951) q[3];
sx q[3];
rz(-0.88328254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6898952) q[0];
sx q[0];
rz(-2.8391333) q[0];
sx q[0];
rz(0.46689335) q[0];
rz(-3.0310071) q[1];
sx q[1];
rz(-0.46375912) q[1];
sx q[1];
rz(2.0707524) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0127481) q[0];
sx q[0];
rz(-1.332988) q[0];
sx q[0];
rz(2.0369916) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7532639) q[2];
sx q[2];
rz(-1.2571722) q[2];
sx q[2];
rz(2.8387866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22510281) q[1];
sx q[1];
rz(-1.9472633) q[1];
sx q[1];
rz(-0.55209362) q[1];
rz(-2.4399906) q[3];
sx q[3];
rz(-2.1013341) q[3];
sx q[3];
rz(-2.4367133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.11833) q[2];
sx q[2];
rz(-1.3619962) q[2];
sx q[2];
rz(-2.2892717) q[2];
rz(0.037633745) q[3];
sx q[3];
rz(-1.2331542) q[3];
sx q[3];
rz(-0.5948624) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1325876) q[0];
sx q[0];
rz(-1.1786893) q[0];
sx q[0];
rz(-1.8527385) q[0];
rz(-1.5124849) q[1];
sx q[1];
rz(-2.0178724) q[1];
sx q[1];
rz(-1.1423133) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2347533) q[0];
sx q[0];
rz(-1.2182353) q[0];
sx q[0];
rz(-1.4666739) q[0];
x q[1];
rz(1.283848) q[2];
sx q[2];
rz(-0.21636886) q[2];
sx q[2];
rz(-1.3138101) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76872494) q[1];
sx q[1];
rz(-1.3228251) q[1];
sx q[1];
rz(-0.3526814) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31123881) q[3];
sx q[3];
rz(-1.4025619) q[3];
sx q[3];
rz(-0.5420891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8196572) q[2];
sx q[2];
rz(-1.2677931) q[2];
sx q[2];
rz(-3.0211871) q[2];
rz(-1.985792) q[3];
sx q[3];
rz(-1.5294231) q[3];
sx q[3];
rz(0.22404484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8026546) q[0];
sx q[0];
rz(-1.3405565) q[0];
sx q[0];
rz(1.4539723) q[0];
rz(-1.5441719) q[1];
sx q[1];
rz(-1.6463966) q[1];
sx q[1];
rz(2.6905751) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2149725) q[0];
sx q[0];
rz(-2.8344013) q[0];
sx q[0];
rz(-0.19812576) q[0];
x q[1];
rz(2.1233929) q[2];
sx q[2];
rz(-2.7180053) q[2];
sx q[2];
rz(-0.46679631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4828795) q[1];
sx q[1];
rz(-0.45156839) q[1];
sx q[1];
rz(0.2803448) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5805832) q[3];
sx q[3];
rz(-2.3811445) q[3];
sx q[3];
rz(1.1475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1958127) q[2];
sx q[2];
rz(-1.4142298) q[2];
sx q[2];
rz(1.7808524) q[2];
rz(-2.1360548) q[3];
sx q[3];
rz(-1.1755627) q[3];
sx q[3];
rz(-1.7520693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239531) q[0];
sx q[0];
rz(-0.12549505) q[0];
sx q[0];
rz(2.4355167) q[0];
rz(-1.5977244) q[1];
sx q[1];
rz(-1.4223301) q[1];
sx q[1];
rz(-2.2854038) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34205758) q[0];
sx q[0];
rz(-1.5634131) q[0];
sx q[0];
rz(-1.5902341) q[0];
rz(0.43120877) q[2];
sx q[2];
rz(-2.2715306) q[2];
sx q[2];
rz(1.7675811) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4256712) q[1];
sx q[1];
rz(-1.7737264) q[1];
sx q[1];
rz(2.5589866) q[1];
rz(3.007902) q[3];
sx q[3];
rz(-1.0858337) q[3];
sx q[3];
rz(0.07459379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2989444) q[2];
sx q[2];
rz(-1.4010669) q[2];
sx q[2];
rz(-2.974158) q[2];
rz(-1.5873448) q[3];
sx q[3];
rz(-0.7730248) q[3];
sx q[3];
rz(-1.0003482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7635968) q[0];
sx q[0];
rz(-1.6908228) q[0];
sx q[0];
rz(-0.3717306) q[0];
rz(-1.4338214) q[1];
sx q[1];
rz(-1.5377518) q[1];
sx q[1];
rz(0.30002123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2903295) q[0];
sx q[0];
rz(-1.4545049) q[0];
sx q[0];
rz(-2.854748) q[0];
rz(-pi) q[1];
rz(-2.5644767) q[2];
sx q[2];
rz(-1.5367998) q[2];
sx q[2];
rz(-2.196072) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.65574291) q[1];
sx q[1];
rz(-2.2169737) q[1];
sx q[1];
rz(1.7263401) q[1];
rz(-pi) q[2];
rz(0.86784243) q[3];
sx q[3];
rz(-1.6658837) q[3];
sx q[3];
rz(-2.8114708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6649449) q[2];
sx q[2];
rz(-2.6740394) q[2];
sx q[2];
rz(-2.7139968) q[2];
rz(1.556373) q[3];
sx q[3];
rz(-2.0245602) q[3];
sx q[3];
rz(-2.3868886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72117358) q[0];
sx q[0];
rz(-0.54062802) q[0];
sx q[0];
rz(-2.6721201) q[0];
rz(0.92533127) q[1];
sx q[1];
rz(-2.053849) q[1];
sx q[1];
rz(-0.023177711) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3566475) q[0];
sx q[0];
rz(-1.0482422) q[0];
sx q[0];
rz(1.3808668) q[0];
rz(2.542664) q[2];
sx q[2];
rz(-0.99236503) q[2];
sx q[2];
rz(-2.3430062) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9564445) q[1];
sx q[1];
rz(-1.7352163) q[1];
sx q[1];
rz(1.3944574) q[1];
x q[2];
rz(-1.3599102) q[3];
sx q[3];
rz(-2.4004705) q[3];
sx q[3];
rz(1.1531342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2716486) q[2];
sx q[2];
rz(-1.2056489) q[2];
sx q[2];
rz(-0.49087697) q[2];
rz(-1.0902181) q[3];
sx q[3];
rz(-0.96929437) q[3];
sx q[3];
rz(-0.31392613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.28521095) q[0];
sx q[0];
rz(-0.87775341) q[0];
sx q[0];
rz(-2.0650771) q[0];
rz(-0.27676997) q[1];
sx q[1];
rz(-0.77817398) q[1];
sx q[1];
rz(-1.4336817) q[1];
rz(0.68885531) q[2];
sx q[2];
rz(-1.8102472) q[2];
sx q[2];
rz(2.3285338) q[2];
rz(-3.1076486) q[3];
sx q[3];
rz(-2.6068519) q[3];
sx q[3];
rz(3.1153771) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
