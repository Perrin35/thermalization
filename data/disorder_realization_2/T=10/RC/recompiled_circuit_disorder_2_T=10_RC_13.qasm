OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4484654) q[0];
sx q[0];
rz(-2.6187596) q[0];
sx q[0];
rz(0.62358207) q[0];
rz(-2.8514255) q[1];
sx q[1];
rz(-0.71915141) q[1];
sx q[1];
rz(-2.6410988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60808027) q[0];
sx q[0];
rz(-2.170616) q[0];
sx q[0];
rz(1.5661201) q[0];
rz(-0.093437336) q[2];
sx q[2];
rz(-1.4215901) q[2];
sx q[2];
rz(1.1364394) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0296368) q[1];
sx q[1];
rz(-1.2458548) q[1];
sx q[1];
rz(-2.3287661) q[1];
x q[2];
rz(-1.1327098) q[3];
sx q[3];
rz(-1.4287018) q[3];
sx q[3];
rz(-0.70664584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7913251) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(-1.2228489) q[2];
rz(1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-0.97035113) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704849) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(2.0626542) q[0];
rz(-1.7547912) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-2.4761377) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42940419) q[0];
sx q[0];
rz(-2.5948338) q[0];
sx q[0];
rz(2.5736546) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7416679) q[2];
sx q[2];
rz(-1.2534281) q[2];
sx q[2];
rz(-0.7764118) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3033777) q[1];
sx q[1];
rz(-1.3347374) q[1];
sx q[1];
rz(-2.6462376) q[1];
x q[2];
rz(-2.2802248) q[3];
sx q[3];
rz(-1.5449761) q[3];
sx q[3];
rz(0.65755075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5370496) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(3.0664505) q[2];
rz(1.4705307) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(-0.69141928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55494088) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(-1.8485908) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(1.741515) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83126691) q[0];
sx q[0];
rz(-1.184549) q[0];
sx q[0];
rz(-1.2826184) q[0];
x q[1];
rz(1.4351575) q[2];
sx q[2];
rz(-1.3856158) q[2];
sx q[2];
rz(0.26091012) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6058265) q[1];
sx q[1];
rz(-1.2488135) q[1];
sx q[1];
rz(0.95226007) q[1];
rz(-pi) q[2];
rz(-1.0812976) q[3];
sx q[3];
rz(-0.34265095) q[3];
sx q[3];
rz(-0.12658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.65163461) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(-1.970132) q[3];
sx q[3];
rz(-1.4843342) q[3];
sx q[3];
rz(-1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3477429) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(1.6963652) q[0];
rz(-1.4472648) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(2.7935374) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845683) q[0];
sx q[0];
rz(-0.94385249) q[0];
sx q[0];
rz(-1.2059962) q[0];
x q[1];
rz(1.1966755) q[2];
sx q[2];
rz(-1.4305563) q[2];
sx q[2];
rz(-1.8082878) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0192249) q[1];
sx q[1];
rz(-0.61673635) q[1];
sx q[1];
rz(1.7255746) q[1];
x q[2];
rz(-0.69339852) q[3];
sx q[3];
rz(-0.20892538) q[3];
sx q[3];
rz(2.8914176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.41670123) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(-0.42281881) q[2];
rz(2.4041798) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(3.056934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(-0.53043956) q[0];
rz(0.92492217) q[1];
sx q[1];
rz(-1.568012) q[1];
sx q[1];
rz(-1.8431429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0153774) q[0];
sx q[0];
rz(-1.4091638) q[0];
sx q[0];
rz(2.9956908) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2735882) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(-2.1538018) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7057719) q[1];
sx q[1];
rz(-0.37322361) q[1];
sx q[1];
rz(-2.865764) q[1];
rz(-pi) q[2];
rz(2.2523746) q[3];
sx q[3];
rz(-2.4272356) q[3];
sx q[3];
rz(1.0970955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9412781) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(-0.47362622) q[2];
rz(3.04223) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6376003) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(-0.9978869) q[0];
rz(-2.2672794) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(2.6748437) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40690639) q[0];
sx q[0];
rz(-2.2285301) q[0];
sx q[0];
rz(1.8440194) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42748638) q[2];
sx q[2];
rz(-1.3689405) q[2];
sx q[2];
rz(1.0496548) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7044428) q[1];
sx q[1];
rz(-1.721742) q[1];
sx q[1];
rz(1.6187861) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.407133) q[3];
sx q[3];
rz(-1.6657077) q[3];
sx q[3];
rz(1.6223736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2816887) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(1.5768645) q[2];
rz(-0.62670296) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(-2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307813) q[0];
sx q[0];
rz(-0.89650506) q[0];
sx q[0];
rz(-2.7217641) q[0];
rz(2.9201674) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(-1.5484757) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683068) q[0];
sx q[0];
rz(-1.5791897) q[0];
sx q[0];
rz(-3.054137) q[0];
x q[1];
rz(2.1599342) q[2];
sx q[2];
rz(-0.33827153) q[2];
sx q[2];
rz(-0.32064082) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6429813) q[1];
sx q[1];
rz(-1.4652518) q[1];
sx q[1];
rz(-1.5555744) q[1];
rz(-2.1920131) q[3];
sx q[3];
rz(-0.55286828) q[3];
sx q[3];
rz(1.0903996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5856813) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(-1.4038203) q[2];
rz(2.3445271) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4306915) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(-2.8523493) q[0];
rz(0.62492433) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(-1.8274868) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24459141) q[0];
sx q[0];
rz(-0.76652157) q[0];
sx q[0];
rz(-0.15890973) q[0];
rz(-pi) q[1];
rz(1.9004702) q[2];
sx q[2];
rz(-0.85061073) q[2];
sx q[2];
rz(2.5230797) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1738759) q[1];
sx q[1];
rz(-1.6248676) q[1];
sx q[1];
rz(1.4500344) q[1];
rz(-pi) q[2];
rz(1.8955599) q[3];
sx q[3];
rz(-2.357558) q[3];
sx q[3];
rz(1.3336381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.73359314) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(-1.4286263) q[2];
rz(-0.96380487) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6190417) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(-1.2458941) q[0];
rz(-3.030581) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(-0.54661173) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7738757) q[0];
sx q[0];
rz(-1.9426632) q[0];
sx q[0];
rz(1.0982151) q[0];
x q[1];
rz(0.39491744) q[2];
sx q[2];
rz(-1.2262218) q[2];
sx q[2];
rz(3.092569) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6706898) q[1];
sx q[1];
rz(-2.7109475) q[1];
sx q[1];
rz(-1.1501269) q[1];
x q[2];
rz(-2.176748) q[3];
sx q[3];
rz(-2.3904739) q[3];
sx q[3];
rz(1.5918819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(1.6938422) q[2];
rz(-2.0041806) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85957134) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(-2.4243673) q[0];
rz(-1.9316797) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(2.4338914) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.087698547) q[0];
sx q[0];
rz(-0.90488926) q[0];
sx q[0];
rz(0.44184394) q[0];
rz(-2.131358) q[2];
sx q[2];
rz(-0.26294225) q[2];
sx q[2];
rz(1.5901142) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0545132) q[1];
sx q[1];
rz(-2.2594249) q[1];
sx q[1];
rz(-2.7282532) q[1];
rz(0.95146146) q[3];
sx q[3];
rz(-1.2776432) q[3];
sx q[3];
rz(-0.64630634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6282965) q[2];
sx q[2];
rz(-2.4847023) q[2];
sx q[2];
rz(2.2383402) q[2];
rz(1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.306504) q[0];
sx q[0];
rz(-2.7764414) q[0];
sx q[0];
rz(2.2055702) q[0];
rz(-2.3256336) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(1.1076526) q[2];
sx q[2];
rz(-1.9911498) q[2];
sx q[2];
rz(-0.1291612) q[2];
rz(3.0872185) q[3];
sx q[3];
rz(-1.5297223) q[3];
sx q[3];
rz(-2.3755304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
