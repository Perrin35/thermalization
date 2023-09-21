OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(-2.5928901) q[0];
sx q[0];
rz(-2.2572416) q[0];
rz(-1.7110775) q[1];
sx q[1];
rz(-0.95354748) q[1];
sx q[1];
rz(-1.5024827) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2384773) q[0];
sx q[0];
rz(-1.7503386) q[0];
sx q[0];
rz(-1.205501) q[0];
rz(-pi) q[1];
rz(-1.5760954) q[2];
sx q[2];
rz(-2.1359348) q[2];
sx q[2];
rz(2.0084755) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.002418) q[1];
sx q[1];
rz(-0.25484172) q[1];
sx q[1];
rz(-2.8789218) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2803335) q[3];
sx q[3];
rz(-1.8464631) q[3];
sx q[3];
rz(1.1112801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(-2.4856429) q[2];
rz(-1.9338699) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2475964) q[0];
sx q[0];
rz(-0.42034724) q[0];
sx q[0];
rz(0.43352747) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(-0.0072335009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9911306) q[0];
sx q[0];
rz(-0.40369243) q[0];
sx q[0];
rz(1.3315721) q[0];
rz(1.8859293) q[2];
sx q[2];
rz(-1.8506146) q[2];
sx q[2];
rz(-1.2535821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1404328) q[1];
sx q[1];
rz(-1.7908485) q[1];
sx q[1];
rz(2.7697255) q[1];
x q[2];
rz(2.6633127) q[3];
sx q[3];
rz(-0.96049958) q[3];
sx q[3];
rz(-1.5680716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6146415) q[2];
sx q[2];
rz(-2.3336637) q[2];
sx q[2];
rz(-2.4439404) q[2];
rz(3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(2.8377623) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4678629) q[0];
sx q[0];
rz(-2.2255852) q[0];
sx q[0];
rz(-1.7720222) q[0];
rz(-1.2415775) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(0.26161584) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25900349) q[0];
sx q[0];
rz(-3.1212174) q[0];
sx q[0];
rz(-2.0632319) q[0];
rz(-pi) q[1];
rz(0.82654731) q[2];
sx q[2];
rz(-1.4790223) q[2];
sx q[2];
rz(-1.9020136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6690327) q[1];
sx q[1];
rz(-1.4497888) q[1];
sx q[1];
rz(-0.78305142) q[1];
x q[2];
rz(1.5393799) q[3];
sx q[3];
rz(-1.0580225) q[3];
sx q[3];
rz(-2.153742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6802784) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(-2.938081) q[2];
rz(0.92173785) q[3];
sx q[3];
rz(-1.2676055) q[3];
sx q[3];
rz(-2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5638567) q[0];
sx q[0];
rz(-1.5699566) q[0];
rz(2.1381901) q[1];
sx q[1];
rz(-1.827821) q[1];
sx q[1];
rz(1.8932231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49444775) q[0];
sx q[0];
rz(-1.6361423) q[0];
sx q[0];
rz(-1.4739743) q[0];
rz(2.4291971) q[2];
sx q[2];
rz(-2.8678896) q[2];
sx q[2];
rz(-1.7143539) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8716988) q[1];
sx q[1];
rz(-1.4596241) q[1];
sx q[1];
rz(-0.65472366) q[1];
x q[2];
rz(1.3094041) q[3];
sx q[3];
rz(-0.86463118) q[3];
sx q[3];
rz(-2.1690926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0288329) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(-1.933243) q[3];
sx q[3];
rz(-1.2669867) q[3];
sx q[3];
rz(-2.3560431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1355302) q[0];
sx q[0];
rz(-2.0879789) q[0];
sx q[0];
rz(0.77520448) q[0];
rz(2.7397621) q[1];
sx q[1];
rz(-2.1907175) q[1];
sx q[1];
rz(-2.2391589) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8309801) q[0];
sx q[0];
rz(-1.6164301) q[0];
sx q[0];
rz(-0.63953416) q[0];
x q[1];
rz(-0.4544223) q[2];
sx q[2];
rz(-1.3497958) q[2];
sx q[2];
rz(-0.99046747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93859766) q[1];
sx q[1];
rz(-1.0293759) q[1];
sx q[1];
rz(2.2375537) q[1];
rz(1.2540713) q[3];
sx q[3];
rz(-1.5734908) q[3];
sx q[3];
rz(0.56841422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.19501413) q[2];
sx q[2];
rz(-0.48823753) q[2];
sx q[2];
rz(-1.1966594) q[2];
rz(1.442391) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8114132) q[0];
sx q[0];
rz(-0.17689642) q[0];
sx q[0];
rz(-0.50317558) q[0];
rz(-1.6852089) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(-2.9398289) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0148894) q[0];
sx q[0];
rz(-3.0198583) q[0];
sx q[0];
rz(-0.99862167) q[0];
rz(2.3927275) q[2];
sx q[2];
rz(-0.61486926) q[2];
sx q[2];
rz(-0.90235898) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0338248) q[1];
sx q[1];
rz(-2.0467842) q[1];
sx q[1];
rz(-2.0460709) q[1];
rz(1.1062578) q[3];
sx q[3];
rz(-0.89336508) q[3];
sx q[3];
rz(0.57074742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9266944) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(-0.8231419) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(2.2657623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.293752) q[0];
sx q[0];
rz(-2.9822615) q[0];
sx q[0];
rz(-3.0840432) q[0];
rz(-1.4808222) q[1];
sx q[1];
rz(-1.3596423) q[1];
sx q[1];
rz(0.94271359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38521117) q[0];
sx q[0];
rz(-0.74059534) q[0];
sx q[0];
rz(-2.5368607) q[0];
rz(-pi) q[1];
rz(-0.16072388) q[2];
sx q[2];
rz(-1.3458999) q[2];
sx q[2];
rz(2.7149534) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47989935) q[1];
sx q[1];
rz(-2.6091895) q[1];
sx q[1];
rz(0.9049306) q[1];
rz(-pi) q[2];
rz(0.25992486) q[3];
sx q[3];
rz(-2.0689031) q[3];
sx q[3];
rz(-0.92344027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1192347) q[2];
sx q[2];
rz(-1.0417577) q[2];
sx q[2];
rz(-0.38267246) q[2];
rz(2.102397) q[3];
sx q[3];
rz(-1.8363876) q[3];
sx q[3];
rz(2.1634845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7169749) q[0];
sx q[0];
rz(-0.096352339) q[0];
sx q[0];
rz(-2.8714645) q[0];
rz(0.62942901) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(-2.8576635) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086872452) q[0];
sx q[0];
rz(-0.76612872) q[0];
sx q[0];
rz(-0.92673577) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4242886) q[2];
sx q[2];
rz(-0.97587817) q[2];
sx q[2];
rz(2.811424) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64360196) q[1];
sx q[1];
rz(-2.0703348) q[1];
sx q[1];
rz(2.1010146) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.038589434) q[3];
sx q[3];
rz(-2.4814197) q[3];
sx q[3];
rz(2.431228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5876864) q[2];
sx q[2];
rz(-2.9710785) q[2];
sx q[2];
rz(-1.930687) q[2];
rz(2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(-2.5134145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07638409) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(-0.22924766) q[0];
rz(0.30300888) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(1.4607666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630336) q[0];
sx q[0];
rz(-1.8543108) q[0];
sx q[0];
rz(-1.5525596) q[0];
rz(-pi) q[1];
x q[1];
rz(2.714459) q[2];
sx q[2];
rz(-2.7286227) q[2];
sx q[2];
rz(2.4954748) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2158828) q[1];
sx q[1];
rz(-2.6162418) q[1];
sx q[1];
rz(-0.044563091) q[1];
rz(2.7542354) q[3];
sx q[3];
rz(-0.82327561) q[3];
sx q[3];
rz(-2.2203317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4298657) q[2];
sx q[2];
rz(-1.9248328) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(2.7097278) q[3];
sx q[3];
rz(-0.46468195) q[3];
sx q[3];
rz(2.4263884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.573695) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(-2.8046872) q[0];
rz(-2.9341872) q[1];
sx q[1];
rz(-2.1131056) q[1];
sx q[1];
rz(-0.38063231) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.80523) q[0];
sx q[0];
rz(-0.21348937) q[0];
sx q[0];
rz(1.9070894) q[0];
rz(-pi) q[1];
rz(2.400488) q[2];
sx q[2];
rz(-1.1198824) q[2];
sx q[2];
rz(-1.2620743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3818647) q[1];
sx q[1];
rz(-1.7654395) q[1];
sx q[1];
rz(-1.0641644) q[1];
rz(-2.5638644) q[3];
sx q[3];
rz(-1.9247705) q[3];
sx q[3];
rz(-0.83565328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(2.005119) q[2];
rz(3.100637) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(-1.827318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8052335) q[0];
sx q[0];
rz(-0.90538607) q[0];
sx q[0];
rz(-0.3120099) q[0];
rz(2.1144755) q[1];
sx q[1];
rz(-1.2925016) q[1];
sx q[1];
rz(2.1137994) q[1];
rz(1.8006051) q[2];
sx q[2];
rz(-1.8125712) q[2];
sx q[2];
rz(-0.3704091) q[2];
rz(-2.9410578) q[3];
sx q[3];
rz(-2.0024828) q[3];
sx q[3];
rz(1.9867292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];