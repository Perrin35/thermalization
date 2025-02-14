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
rz(0.62975878) q[0];
sx q[0];
rz(3.44343) q[0];
sx q[0];
rz(11.259196) q[0];
rz(2.6131926) q[1];
sx q[1];
rz(-2.3101248) q[1];
sx q[1];
rz(2.6822579) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0889657) q[0];
sx q[0];
rz(-1.8267434) q[0];
sx q[0];
rz(0.050764485) q[0];
rz(-pi) q[1];
rz(1.2940965) q[2];
sx q[2];
rz(-1.7613698) q[2];
sx q[2];
rz(2.0399567) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3078008) q[1];
sx q[1];
rz(-1.3978617) q[1];
sx q[1];
rz(0.61638919) q[1];
rz(-2.2396022) q[3];
sx q[3];
rz(-1.588527) q[3];
sx q[3];
rz(-0.56541857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1186195) q[2];
sx q[2];
rz(-1.8034673) q[2];
sx q[2];
rz(2.5845134) q[2];
rz(2.6856375) q[3];
sx q[3];
rz(-2.3796701) q[3];
sx q[3];
rz(2.5383811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.93422455) q[0];
sx q[0];
rz(-2.4226483) q[0];
sx q[0];
rz(1.3907322) q[0];
rz(-1.8234183) q[1];
sx q[1];
rz(-1.428182) q[1];
sx q[1];
rz(2.9255829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073450449) q[0];
sx q[0];
rz(-1.2552069) q[0];
sx q[0];
rz(0.95629779) q[0];
rz(0.26311263) q[2];
sx q[2];
rz(-2.7091658) q[2];
sx q[2];
rz(-0.95460923) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29874814) q[1];
sx q[1];
rz(-2.2262888) q[1];
sx q[1];
rz(2.4821494) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9253821) q[3];
sx q[3];
rz(-0.56963339) q[3];
sx q[3];
rz(-2.5550752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7184489) q[2];
sx q[2];
rz(-0.37280145) q[2];
sx q[2];
rz(-0.049840363) q[2];
rz(0.54346624) q[3];
sx q[3];
rz(-1.1246357) q[3];
sx q[3];
rz(2.0487093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8922888) q[0];
sx q[0];
rz(-2.0396621) q[0];
sx q[0];
rz(-0.9683384) q[0];
rz(0.020542055) q[1];
sx q[1];
rz(-1.7612709) q[1];
sx q[1];
rz(-1.4302018) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50043369) q[0];
sx q[0];
rz(-2.4999188) q[0];
sx q[0];
rz(-3.1042751) q[0];
x q[1];
rz(-0.5495785) q[2];
sx q[2];
rz(-2.5039154) q[2];
sx q[2];
rz(0.23070947) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80950981) q[1];
sx q[1];
rz(-0.44594279) q[1];
sx q[1];
rz(-2.6757702) q[1];
rz(-pi) q[2];
rz(2.7647721) q[3];
sx q[3];
rz(-1.6018036) q[3];
sx q[3];
rz(-3.0135807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.7057544) q[2];
sx q[2];
rz(-0.64283723) q[2];
sx q[2];
rz(1.359681) q[2];
rz(2.6333574) q[3];
sx q[3];
rz(-1.619092) q[3];
sx q[3];
rz(-1.9967509) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0559167) q[0];
sx q[0];
rz(-2.8450232) q[0];
sx q[0];
rz(-2.1348409) q[0];
rz(0.83167568) q[1];
sx q[1];
rz(-0.74200231) q[1];
sx q[1];
rz(1.1078018) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9079893) q[0];
sx q[0];
rz(-2.395326) q[0];
sx q[0];
rz(2.1265592) q[0];
x q[1];
rz(1.8227656) q[2];
sx q[2];
rz(-1.2059847) q[2];
sx q[2];
rz(3.0598726) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7101871) q[1];
sx q[1];
rz(-1.5509938) q[1];
sx q[1];
rz(1.8777678) q[1];
rz(1.2899701) q[3];
sx q[3];
rz(-2.0005595) q[3];
sx q[3];
rz(-2.6971779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45336777) q[2];
sx q[2];
rz(-1.0658762) q[2];
sx q[2];
rz(-0.31309703) q[2];
rz(-2.744216) q[3];
sx q[3];
rz(-1.4704967) q[3];
sx q[3];
rz(-0.1853005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.6322286) q[0];
sx q[0];
rz(-2.2780184) q[0];
sx q[0];
rz(2.2160227) q[0];
rz(-0.46982345) q[1];
sx q[1];
rz(-1.8269822) q[1];
sx q[1];
rz(2.125461) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52262703) q[0];
sx q[0];
rz(-0.62331182) q[0];
sx q[0];
rz(2.3254184) q[0];
rz(-pi) q[1];
rz(-2.2158898) q[2];
sx q[2];
rz(-0.4302372) q[2];
sx q[2];
rz(1.0412585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4279856) q[1];
sx q[1];
rz(-1.3283037) q[1];
sx q[1];
rz(-2.9338783) q[1];
rz(3.1096439) q[3];
sx q[3];
rz(-1.2476139) q[3];
sx q[3];
rz(1.1695216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.65285811) q[2];
sx q[2];
rz(-2.8159339) q[2];
sx q[2];
rz(2.6540836) q[2];
rz(-0.89753914) q[3];
sx q[3];
rz(-1.2582015) q[3];
sx q[3];
rz(-1.0645617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2010736) q[0];
sx q[0];
rz(-2.2036393) q[0];
sx q[0];
rz(0.67365375) q[0];
rz(-1.9081839) q[1];
sx q[1];
rz(-1.0181095) q[1];
sx q[1];
rz(-0.76603755) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8476794) q[0];
sx q[0];
rz(-1.3013757) q[0];
sx q[0];
rz(1.040017) q[0];
x q[1];
rz(-0.38389194) q[2];
sx q[2];
rz(-2.9523179) q[2];
sx q[2];
rz(-1.4658907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1704857) q[1];
sx q[1];
rz(-2.0506713) q[1];
sx q[1];
rz(2.9779153) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4333731) q[3];
sx q[3];
rz(-1.2677464) q[3];
sx q[3];
rz(2.0772484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.96887863) q[2];
sx q[2];
rz(-1.4977027) q[2];
sx q[2];
rz(2.3806351) q[2];
rz(-0.40192762) q[3];
sx q[3];
rz(-2.8765078) q[3];
sx q[3];
rz(-2.5529805) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48259476) q[0];
sx q[0];
rz(-0.75356475) q[0];
sx q[0];
rz(1.6931417) q[0];
rz(0.50085577) q[1];
sx q[1];
rz(-1.8468937) q[1];
sx q[1];
rz(-0.81333152) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9324009) q[0];
sx q[0];
rz(-1.4228357) q[0];
sx q[0];
rz(2.7514821) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3766889) q[2];
sx q[2];
rz(-0.55667215) q[2];
sx q[2];
rz(2.106032) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2709368) q[1];
sx q[1];
rz(-1.7084048) q[1];
sx q[1];
rz(0.44160053) q[1];
x q[2];
rz(-1.6930831) q[3];
sx q[3];
rz(-1.5124195) q[3];
sx q[3];
rz(-0.98298847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8616051) q[2];
sx q[2];
rz(-2.5364752) q[2];
sx q[2];
rz(2.81847) q[2];
rz(1.4019639) q[3];
sx q[3];
rz(-1.2819382) q[3];
sx q[3];
rz(-1.3824979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0793656) q[0];
sx q[0];
rz(-1.0492188) q[0];
sx q[0];
rz(0.76612377) q[0];
rz(-0.8194204) q[1];
sx q[1];
rz(-2.1111646) q[1];
sx q[1];
rz(1.6105509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4774276) q[0];
sx q[0];
rz(-2.9026845) q[0];
sx q[0];
rz(-2.8498123) q[0];
x q[1];
rz(1.3513597) q[2];
sx q[2];
rz(-0.92715741) q[2];
sx q[2];
rz(2.6257169) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.69680113) q[1];
sx q[1];
rz(-2.2292622) q[1];
sx q[1];
rz(-0.80689238) q[1];
x q[2];
rz(0.29115486) q[3];
sx q[3];
rz(-1.8829573) q[3];
sx q[3];
rz(-0.96033421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7649585) q[2];
sx q[2];
rz(-0.19516334) q[2];
sx q[2];
rz(-2.7435319) q[2];
rz(0.6063439) q[3];
sx q[3];
rz(-1.5747993) q[3];
sx q[3];
rz(-0.91355598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23676087) q[0];
sx q[0];
rz(-0.30192152) q[0];
sx q[0];
rz(-0.52445573) q[0];
rz(2.4728788) q[1];
sx q[1];
rz(-1.5293744) q[1];
sx q[1];
rz(1.761577) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83630122) q[0];
sx q[0];
rz(-0.99812767) q[0];
sx q[0];
rz(0.051866626) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0933268) q[2];
sx q[2];
rz(-1.2984972) q[2];
sx q[2];
rz(1.7881025) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3375597) q[1];
sx q[1];
rz(-1.1566678) q[1];
sx q[1];
rz(-1.1131338) q[1];
x q[2];
rz(-0.53219019) q[3];
sx q[3];
rz(-1.8992153) q[3];
sx q[3];
rz(1.6468954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9809208) q[2];
sx q[2];
rz(-0.74690861) q[2];
sx q[2];
rz(-0.46372947) q[2];
rz(1.7175698) q[3];
sx q[3];
rz(-1.0878891) q[3];
sx q[3];
rz(-3.1080642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.5793107) q[0];
sx q[0];
rz(-0.71826851) q[0];
sx q[0];
rz(2.723519) q[0];
rz(-2.383291) q[1];
sx q[1];
rz(-2.1577991) q[1];
sx q[1];
rz(-1.8565348) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4649974) q[0];
sx q[0];
rz(-2.2187859) q[0];
sx q[0];
rz(1.4761488) q[0];
x q[1];
rz(2.0359155) q[2];
sx q[2];
rz(-2.4659202) q[2];
sx q[2];
rz(-0.0059520324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.549379) q[1];
sx q[1];
rz(-0.58281189) q[1];
sx q[1];
rz(-2.7874376) q[1];
rz(-pi) q[2];
rz(0.22104903) q[3];
sx q[3];
rz(-1.9004824) q[3];
sx q[3];
rz(0.21176618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.16984223) q[2];
sx q[2];
rz(-2.1368133) q[2];
sx q[2];
rz(-0.24610914) q[2];
rz(1.2717815) q[3];
sx q[3];
rz(-1.5747986) q[3];
sx q[3];
rz(1.3625712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1590189) q[0];
sx q[0];
rz(-1.6652501) q[0];
sx q[0];
rz(-0.2022947) q[0];
rz(2.5166439) q[1];
sx q[1];
rz(-2.2849871) q[1];
sx q[1];
rz(-2.7402592) q[1];
rz(-1.779378) q[2];
sx q[2];
rz(-1.7558395) q[2];
sx q[2];
rz(0.21370733) q[2];
rz(2.940964) q[3];
sx q[3];
rz(-1.8798141) q[3];
sx q[3];
rz(0.78165913) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
