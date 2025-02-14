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
rz(-0.023802726) q[0];
sx q[0];
rz(-2.0355712) q[0];
sx q[0];
rz(0.7769146) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(2.1652752) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1084749) q[0];
sx q[0];
rz(-0.43370789) q[0];
sx q[0];
rz(1.2483622) q[0];
rz(2.5778779) q[2];
sx q[2];
rz(-1.6275121) q[2];
sx q[2];
rz(2.4617755) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0529671) q[1];
sx q[1];
rz(-1.2994811) q[1];
sx q[1];
rz(2.0120828) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5209274) q[3];
sx q[3];
rz(-2.5679776) q[3];
sx q[3];
rz(-2.916934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.50600791) q[2];
sx q[2];
rz(-0.053746544) q[2];
sx q[2];
rz(0.40679833) q[2];
rz(2.9721416) q[3];
sx q[3];
rz(-0.5295161) q[3];
sx q[3];
rz(-1.0725526) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605543) q[0];
sx q[0];
rz(-2.9091703) q[0];
sx q[0];
rz(0.0037923092) q[0];
rz(-3.0637528) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(-2.8357764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90214848) q[0];
sx q[0];
rz(-1.7551219) q[0];
sx q[0];
rz(-1.0499766) q[0];
x q[1];
rz(-2.9980837) q[2];
sx q[2];
rz(-2.3845551) q[2];
sx q[2];
rz(-0.17998634) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0111573) q[1];
sx q[1];
rz(-0.88825916) q[1];
sx q[1];
rz(0.52027563) q[1];
rz(-pi) q[2];
rz(0.81021328) q[3];
sx q[3];
rz(-3.0533724) q[3];
sx q[3];
rz(0.42217964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99016142) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(-2.4988417) q[2];
rz(0.07240545) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(-1.7806627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088242315) q[0];
sx q[0];
rz(-0.14277661) q[0];
sx q[0];
rz(-2.9852168) q[0];
rz(-0.0414255) q[1];
sx q[1];
rz(-0.62774575) q[1];
sx q[1];
rz(1.590439) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8122221) q[0];
sx q[0];
rz(-0.58766627) q[0];
sx q[0];
rz(-0.37826041) q[0];
x q[1];
rz(1.5029491) q[2];
sx q[2];
rz(-0.91582752) q[2];
sx q[2];
rz(-0.52815765) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2385905) q[1];
sx q[1];
rz(-0.5410453) q[1];
sx q[1];
rz(0.21958406) q[1];
rz(2.2699192) q[3];
sx q[3];
rz(-2.1101885) q[3];
sx q[3];
rz(-1.9339069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40858832) q[2];
sx q[2];
rz(-2.6027347) q[2];
sx q[2];
rz(-0.020922529) q[2];
rz(2.9544592) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(0.11370295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1784096) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(2.7030429) q[0];
rz(1.6167538) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(0.24510342) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.995718) q[0];
sx q[0];
rz(-0.85332131) q[0];
sx q[0];
rz(0.11086734) q[0];
rz(-pi) q[1];
rz(-2.7510277) q[2];
sx q[2];
rz(-2.7465944) q[2];
sx q[2];
rz(2.4633138) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72232258) q[1];
sx q[1];
rz(-2.5825204) q[1];
sx q[1];
rz(-1.30617) q[1];
rz(-pi) q[2];
rz(1.887759) q[3];
sx q[3];
rz(-0.84268236) q[3];
sx q[3];
rz(-1.9469946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(-2.8098246) q[2];
rz(0.48745421) q[3];
sx q[3];
rz(-1.0468227) q[3];
sx q[3];
rz(-2.148518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29193923) q[0];
sx q[0];
rz(-1.6878457) q[0];
sx q[0];
rz(-0.77350235) q[0];
rz(-1.9603112) q[1];
sx q[1];
rz(-0.14207323) q[1];
sx q[1];
rz(1.7519417) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6365125) q[0];
sx q[0];
rz(-0.55547041) q[0];
sx q[0];
rz(2.4019) q[0];
x q[1];
rz(0.85948617) q[2];
sx q[2];
rz(-2.0469249) q[2];
sx q[2];
rz(-1.7993594) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.90926814) q[1];
sx q[1];
rz(-1.0977912) q[1];
sx q[1];
rz(-1.478341) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7520653) q[3];
sx q[3];
rz(-1.3366404) q[3];
sx q[3];
rz(2.2088946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2191849) q[2];
sx q[2];
rz(-1.2097404) q[2];
sx q[2];
rz(-0.47214559) q[2];
rz(1.2989429) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(-0.81056547) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2209114) q[0];
sx q[0];
rz(-0.26316106) q[0];
sx q[0];
rz(0.26350185) q[0];
rz(1.1031411) q[1];
sx q[1];
rz(-1.3145072) q[1];
sx q[1];
rz(-2.7679494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6275245) q[0];
sx q[0];
rz(-1.4980982) q[0];
sx q[0];
rz(0.060427314) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0194026) q[2];
sx q[2];
rz(-2.2940679) q[2];
sx q[2];
rz(-0.10553372) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9615104) q[1];
sx q[1];
rz(-1.0133378) q[1];
sx q[1];
rz(0.37661676) q[1];
rz(-pi) q[2];
rz(2.5534036) q[3];
sx q[3];
rz(-1.08687) q[3];
sx q[3];
rz(-1.5744792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0282447) q[2];
sx q[2];
rz(-0.17013203) q[2];
sx q[2];
rz(-2.587758) q[2];
rz(-1.3977741) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(-0.3127313) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58436191) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(1.9118017) q[0];
rz(2.9025485) q[1];
sx q[1];
rz(-1.5115279) q[1];
sx q[1];
rz(-0.30034932) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88278937) q[0];
sx q[0];
rz(-2.9414294) q[0];
sx q[0];
rz(-0.61997719) q[0];
rz(-pi) q[1];
rz(2.5744252) q[2];
sx q[2];
rz(-2.3544899) q[2];
sx q[2];
rz(1.2445104) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7869889) q[1];
sx q[1];
rz(-1.5704078) q[1];
sx q[1];
rz(-1.555519) q[1];
rz(-pi) q[2];
rz(3.1170397) q[3];
sx q[3];
rz(-0.85806393) q[3];
sx q[3];
rz(-2.7803382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0099237) q[2];
sx q[2];
rz(-1.4574304) q[2];
sx q[2];
rz(0.24492502) q[2];
rz(0.51982546) q[3];
sx q[3];
rz(-2.2782785) q[3];
sx q[3];
rz(-0.68827099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7898665) q[0];
sx q[0];
rz(-2.7930197) q[0];
sx q[0];
rz(-1.1630195) q[0];
rz(0.066896833) q[1];
sx q[1];
rz(-1.4935378) q[1];
sx q[1];
rz(-2.1287207) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.430535) q[0];
sx q[0];
rz(-1.0037046) q[0];
sx q[0];
rz(2.5398769) q[0];
x q[1];
rz(-1.2946473) q[2];
sx q[2];
rz(-2.5564402) q[2];
sx q[2];
rz(0.25979751) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2248963) q[1];
sx q[1];
rz(-0.46474296) q[1];
sx q[1];
rz(2.3791989) q[1];
rz(-0.68655218) q[3];
sx q[3];
rz(-1.2408537) q[3];
sx q[3];
rz(0.8984962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.11671994) q[2];
sx q[2];
rz(-2.1634384) q[2];
sx q[2];
rz(-0.32279521) q[2];
rz(-0.60574496) q[3];
sx q[3];
rz(-0.79234684) q[3];
sx q[3];
rz(-0.34887031) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054319687) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(-0.012454575) q[0];
rz(-0.74673486) q[1];
sx q[1];
rz(-2.2181999) q[1];
sx q[1];
rz(-2.8616203) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075172193) q[0];
sx q[0];
rz(-3.018258) q[0];
sx q[0];
rz(1.146011) q[0];
rz(-pi) q[1];
rz(-0.32692636) q[2];
sx q[2];
rz(-1.9844311) q[2];
sx q[2];
rz(-1.3437611) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2586883) q[1];
sx q[1];
rz(-1.8778442) q[1];
sx q[1];
rz(-0.76489246) q[1];
x q[2];
rz(2.1567508) q[3];
sx q[3];
rz(-0.42736125) q[3];
sx q[3];
rz(0.40287429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.81194699) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(-2.9296056) q[2];
rz(0.82344615) q[3];
sx q[3];
rz(-1.6971089) q[3];
sx q[3];
rz(-2.8823891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9987746) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(2.4488191) q[0];
rz(0.57299262) q[1];
sx q[1];
rz(-1.3447821) q[1];
sx q[1];
rz(-0.43100345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8739024) q[0];
sx q[0];
rz(-2.5866716) q[0];
sx q[0];
rz(-0.11124994) q[0];
rz(-pi) q[1];
rz(-2.0153322) q[2];
sx q[2];
rz(-1.017414) q[2];
sx q[2];
rz(2.0972507) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5332408) q[1];
sx q[1];
rz(-2.2679288) q[1];
sx q[1];
rz(1.240154) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0968571) q[3];
sx q[3];
rz(-2.5383679) q[3];
sx q[3];
rz(-0.67765498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4219605) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(2.5813622) q[2];
rz(2.6719921) q[3];
sx q[3];
rz(-2.738939) q[3];
sx q[3];
rz(2.4277021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2847168) q[0];
sx q[0];
rz(-1.4061883) q[0];
sx q[0];
rz(-1.0549369) q[0];
rz(0.84125413) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(2.8259059) q[2];
sx q[2];
rz(-1.9236947) q[2];
sx q[2];
rz(-2.535939) q[2];
rz(1.5359405) q[3];
sx q[3];
rz(-2.4823454) q[3];
sx q[3];
rz(-2.1128826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
