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
rz(0.99211168) q[0];
sx q[0];
rz(3.8881128) q[0];
sx q[0];
rz(9.9915656) q[0];
rz(1.2504638) q[1];
sx q[1];
rz(-1.992978) q[1];
sx q[1];
rz(2.221938) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56432589) q[0];
sx q[0];
rz(-2.3503605) q[0];
sx q[0];
rz(-1.7971695) q[0];
rz(0.72534277) q[2];
sx q[2];
rz(-1.9646461) q[2];
sx q[2];
rz(-2.3462636) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2055605) q[1];
sx q[1];
rz(-1.524393) q[1];
sx q[1];
rz(-0.64179365) q[1];
rz(0.23969526) q[3];
sx q[3];
rz(-2.5928232) q[3];
sx q[3];
rz(-1.5645998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93996843) q[2];
sx q[2];
rz(-2.2087966) q[2];
sx q[2];
rz(2.0882108) q[2];
rz(0.71422226) q[3];
sx q[3];
rz(-0.37319365) q[3];
sx q[3];
rz(-1.8478954) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6222222) q[0];
sx q[0];
rz(-2.433233) q[0];
sx q[0];
rz(2.569662) q[0];
rz(-1.8427303) q[1];
sx q[1];
rz(-2.3426901) q[1];
sx q[1];
rz(-0.78972185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2970726) q[0];
sx q[0];
rz(-1.6668586) q[0];
sx q[0];
rz(1.2681095) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4873494) q[2];
sx q[2];
rz(-1.6983319) q[2];
sx q[2];
rz(2.5728284) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.29282031) q[1];
sx q[1];
rz(-2.1629984) q[1];
sx q[1];
rz(-1.2315537) q[1];
rz(-pi) q[2];
rz(-0.67482194) q[3];
sx q[3];
rz(-1.1742419) q[3];
sx q[3];
rz(-3.0930072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.85084891) q[2];
sx q[2];
rz(-0.76883832) q[2];
sx q[2];
rz(-1.3624462) q[2];
rz(-0.29081523) q[3];
sx q[3];
rz(-1.7751866) q[3];
sx q[3];
rz(-0.10663685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13880759) q[0];
sx q[0];
rz(-1.0951575) q[0];
sx q[0];
rz(-0.7005257) q[0];
rz(-0.48577148) q[1];
sx q[1];
rz(-2.3120717) q[1];
sx q[1];
rz(-0.22451678) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4493664) q[0];
sx q[0];
rz(-0.3825408) q[0];
sx q[0];
rz(-0.40479779) q[0];
rz(1.4204558) q[2];
sx q[2];
rz(-2.5579961) q[2];
sx q[2];
rz(-0.70202578) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77188801) q[1];
sx q[1];
rz(-1.4955031) q[1];
sx q[1];
rz(0.98702191) q[1];
rz(0.84214476) q[3];
sx q[3];
rz(-1.6557459) q[3];
sx q[3];
rz(-0.85499518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38480467) q[2];
sx q[2];
rz(-0.89185682) q[2];
sx q[2];
rz(3.0925114) q[2];
rz(-0.067042025) q[3];
sx q[3];
rz(-1.9837572) q[3];
sx q[3];
rz(2.4750347) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1753569) q[0];
sx q[0];
rz(-2.5844564) q[0];
sx q[0];
rz(0.019158451) q[0];
rz(-1.7823904) q[1];
sx q[1];
rz(-1.4298341) q[1];
sx q[1];
rz(3.137099) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2990103) q[0];
sx q[0];
rz(-2.0662464) q[0];
sx q[0];
rz(0.61758496) q[0];
x q[1];
rz(2.6439444) q[2];
sx q[2];
rz(-1.9714173) q[2];
sx q[2];
rz(1.7965574) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.055306074) q[1];
sx q[1];
rz(-1.9290566) q[1];
sx q[1];
rz(0.18017027) q[1];
rz(2.9380083) q[3];
sx q[3];
rz(-1.4798375) q[3];
sx q[3];
rz(1.1435777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4517453) q[2];
sx q[2];
rz(-2.0514026) q[2];
sx q[2];
rz(1.6881662) q[2];
rz(1.8870185) q[3];
sx q[3];
rz(-1.5213608) q[3];
sx q[3];
rz(-0.5622676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8101863) q[0];
sx q[0];
rz(-1.2368546) q[0];
sx q[0];
rz(-1.9796665) q[0];
rz(0.76238531) q[1];
sx q[1];
rz(-0.98680174) q[1];
sx q[1];
rz(2.7120178) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7976101) q[0];
sx q[0];
rz(-1.5344041) q[0];
sx q[0];
rz(-1.9790566) q[0];
rz(-pi) q[1];
rz(-2.3829382) q[2];
sx q[2];
rz(-1.3765556) q[2];
sx q[2];
rz(1.615103) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.67744327) q[1];
sx q[1];
rz(-0.92532571) q[1];
sx q[1];
rz(0.42663017) q[1];
rz(-2.2240673) q[3];
sx q[3];
rz(-1.297373) q[3];
sx q[3];
rz(1.3454252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.21542159) q[2];
sx q[2];
rz(-1.2913707) q[2];
sx q[2];
rz(-1.8298836) q[2];
rz(-2.6808776) q[3];
sx q[3];
rz(-2.1381133) q[3];
sx q[3];
rz(0.13319143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.8120414) q[0];
sx q[0];
rz(-2.9819745) q[0];
sx q[0];
rz(2.0352236) q[0];
rz(-0.1144935) q[1];
sx q[1];
rz(-1.0527) q[1];
sx q[1];
rz(0.053650275) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.69218) q[0];
sx q[0];
rz(-0.25568889) q[0];
sx q[0];
rz(1.5935672) q[0];
rz(-pi) q[1];
rz(2.3340763) q[2];
sx q[2];
rz(-1.4984594) q[2];
sx q[2];
rz(1.4405719) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.8369537) q[1];
sx q[1];
rz(-1.0598247) q[1];
sx q[1];
rz(3.0797019) q[1];
rz(2.1838004) q[3];
sx q[3];
rz(-1.0292064) q[3];
sx q[3];
rz(-1.047618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8657118) q[2];
sx q[2];
rz(-1.2578473) q[2];
sx q[2];
rz(-2.6344521) q[2];
rz(1.3745314) q[3];
sx q[3];
rz(-1.5406698) q[3];
sx q[3];
rz(1.1486294) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5116665) q[0];
sx q[0];
rz(-1.1277132) q[0];
sx q[0];
rz(3.1003057) q[0];
rz(2.4033026) q[1];
sx q[1];
rz(-1.9169044) q[1];
sx q[1];
rz(-2.1521177) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51136298) q[0];
sx q[0];
rz(-2.0718899) q[0];
sx q[0];
rz(-0.54176919) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.760354) q[2];
sx q[2];
rz(-1.9206534) q[2];
sx q[2];
rz(0.39385349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.38997) q[1];
sx q[1];
rz(-0.97725673) q[1];
sx q[1];
rz(-0.71243993) q[1];
rz(2.2253898) q[3];
sx q[3];
rz(-1.8713017) q[3];
sx q[3];
rz(2.253807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5993293) q[2];
sx q[2];
rz(-2.0928536) q[2];
sx q[2];
rz(2.2824724) q[2];
rz(3.1298992) q[3];
sx q[3];
rz(-1.499736) q[3];
sx q[3];
rz(-0.11277994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.3149253) q[0];
sx q[0];
rz(-2.5496917) q[0];
sx q[0];
rz(2.9097606) q[0];
rz(2.5595317) q[1];
sx q[1];
rz(-1.8128017) q[1];
sx q[1];
rz(-1.2190762) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7340379) q[0];
sx q[0];
rz(-1.014632) q[0];
sx q[0];
rz(1.708605) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35301669) q[2];
sx q[2];
rz(-1.5215877) q[2];
sx q[2];
rz(-2.4491058) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6775359) q[1];
sx q[1];
rz(-1.4194064) q[1];
sx q[1];
rz(-1.7580126) q[1];
x q[2];
rz(2.3858504) q[3];
sx q[3];
rz(-1.1583987) q[3];
sx q[3];
rz(-2.2568767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1845188) q[2];
sx q[2];
rz(-1.3513214) q[2];
sx q[2];
rz(0.78682023) q[2];
rz(-1.6400853) q[3];
sx q[3];
rz(-0.4314751) q[3];
sx q[3];
rz(1.7155581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1805434) q[0];
sx q[0];
rz(-1.7683872) q[0];
sx q[0];
rz(2.2013262) q[0];
rz(2.6967948) q[1];
sx q[1];
rz(-2.2856789) q[1];
sx q[1];
rz(-0.14642265) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43035591) q[0];
sx q[0];
rz(-2.8161263) q[0];
sx q[0];
rz(-0.40367608) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8249885) q[2];
sx q[2];
rz(-1.3066402) q[2];
sx q[2];
rz(-3.11655) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.47825731) q[1];
sx q[1];
rz(-0.90090776) q[1];
sx q[1];
rz(1.2557593) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10036631) q[3];
sx q[3];
rz(-2.3477738) q[3];
sx q[3];
rz(0.19089261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.08708295) q[2];
sx q[2];
rz(-1.4887709) q[2];
sx q[2];
rz(2.3040237) q[2];
rz(0.93291035) q[3];
sx q[3];
rz(-2.1541903) q[3];
sx q[3];
rz(-2.3605409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26302108) q[0];
sx q[0];
rz(-1.2282547) q[0];
sx q[0];
rz(1.3680869) q[0];
rz(3.0655762) q[1];
sx q[1];
rz(-2.5612505) q[1];
sx q[1];
rz(-2.0215633) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0033542643) q[0];
sx q[0];
rz(-1.8347667) q[0];
sx q[0];
rz(-1.8916983) q[0];
x q[1];
rz(1.2659968) q[2];
sx q[2];
rz(-2.2561361) q[2];
sx q[2];
rz(3.1076269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.489913) q[1];
sx q[1];
rz(-1.8615684) q[1];
sx q[1];
rz(-3.0417697) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.038123) q[3];
sx q[3];
rz(-1.7333687) q[3];
sx q[3];
rz(-1.181447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7398305) q[2];
sx q[2];
rz(-1.4515452) q[2];
sx q[2];
rz(2.8945727) q[2];
rz(-0.57010993) q[3];
sx q[3];
rz(-2.6542122) q[3];
sx q[3];
rz(1.1183687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0982672) q[0];
sx q[0];
rz(-1.7541616) q[0];
sx q[0];
rz(2.7761205) q[0];
rz(3.0430766) q[1];
sx q[1];
rz(-2.5010074) q[1];
sx q[1];
rz(0.88190257) q[1];
rz(0.013507387) q[2];
sx q[2];
rz(-1.1579222) q[2];
sx q[2];
rz(2.9029878) q[2];
rz(2.7612272) q[3];
sx q[3];
rz(-1.3283397) q[3];
sx q[3];
rz(2.0793177) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
