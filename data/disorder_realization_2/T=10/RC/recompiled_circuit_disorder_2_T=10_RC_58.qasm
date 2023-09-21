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
rz(-2.5180106) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(-0.50049385) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96007632) q[0];
sx q[0];
rz(-1.5669364) q[0];
sx q[0];
rz(0.5998248) q[0];
rz(-1.0153158) q[2];
sx q[2];
rz(-0.17586389) q[2];
sx q[2];
rz(1.6989087) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7841255) q[1];
sx q[1];
rz(-0.81175121) q[1];
sx q[1];
rz(2.0264506) q[1];
x q[2];
rz(-2.9849103) q[3];
sx q[3];
rz(-2.0041668) q[3];
sx q[3];
rz(0.79790786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7913251) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(-1.9187437) q[2];
rz(-1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7704849) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(1.0789385) q[0];
rz(-1.7547912) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-2.4761377) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7121885) q[0];
sx q[0];
rz(-0.5467589) q[0];
sx q[0];
rz(2.5736546) q[0];
x q[1];
rz(-0.47768728) q[2];
sx q[2];
rz(-2.7825232) q[2];
sx q[2];
rz(1.2815086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0004955) q[1];
sx q[1];
rz(-2.5971203) q[1];
sx q[1];
rz(-2.673124) q[1];
rz(-pi) q[2];
rz(0.034025107) q[3];
sx q[3];
rz(-0.861654) q[3];
sx q[3];
rz(2.2505086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5370496) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(1.4705307) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55494088) q[0];
sx q[0];
rz(-1.8692769) q[0];
sx q[0];
rz(2.3828322) q[0];
rz(-1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(1.4000777) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5132719) q[0];
sx q[0];
rz(-1.3043881) q[0];
sx q[0];
rz(2.7404286) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6252458) q[2];
sx q[2];
rz(-2.9125104) q[2];
sx q[2];
rz(-2.2428227) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6058265) q[1];
sx q[1];
rz(-1.2488135) q[1];
sx q[1];
rz(2.1893326) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9754144) q[3];
sx q[3];
rz(-1.8718534) q[3];
sx q[3];
rz(-0.64134669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65163461) q[2];
sx q[2];
rz(-0.48214665) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.79384971) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(-1.6963652) q[0];
rz(-1.4472648) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(-2.7935374) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8349583) q[0];
sx q[0];
rz(-1.8638896) q[0];
sx q[0];
rz(-2.4819863) q[0];
rz(2.9910827) q[2];
sx q[2];
rz(-1.2005271) q[2];
sx q[2];
rz(-2.9589257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0192249) q[1];
sx q[1];
rz(-2.5248563) q[1];
sx q[1];
rz(-1.4160181) q[1];
rz(1.4361037) q[3];
sx q[3];
rz(-1.4105984) q[3];
sx q[3];
rz(-0.45405162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7248914) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(-2.7187738) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(-0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2622862) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(-0.53043956) q[0];
rz(-0.92492217) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(1.2984498) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86622483) q[0];
sx q[0];
rz(-0.21731649) q[0];
sx q[0];
rz(-0.84262459) q[0];
rz(-pi) q[1];
rz(1.8680044) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(0.98779087) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4106771) q[1];
sx q[1];
rz(-1.2123322) q[1];
sx q[1];
rz(-1.4645542) q[1];
x q[2];
rz(-2.1634444) q[3];
sx q[3];
rz(-1.9962629) q[3];
sx q[3];
rz(3.0654207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(0.47362622) q[2];
rz(0.099362699) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(0.84053269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50399238) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(2.1437058) q[0];
rz(-2.2672794) q[1];
sx q[1];
rz(-1.0214146) q[1];
sx q[1];
rz(-0.46674892) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3335553) q[0];
sx q[0];
rz(-1.3555962) q[0];
sx q[0];
rz(-2.4654885) q[0];
rz(-2.6830707) q[2];
sx q[2];
rz(-0.47007559) q[2];
sx q[2];
rz(-2.2058861) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14086831) q[1];
sx q[1];
rz(-1.6182401) q[1];
sx q[1];
rz(2.9904757) q[1];
rz(1.407133) q[3];
sx q[3];
rz(-1.4758849) q[3];
sx q[3];
rz(-1.519219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85990396) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(-1.5768645) q[2];
rz(2.5148897) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(-0.68157649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307813) q[0];
sx q[0];
rz(-0.89650506) q[0];
sx q[0];
rz(-2.7217641) q[0];
rz(-0.22142521) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(-1.5484757) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2433462) q[0];
sx q[0];
rz(-1.6582489) q[0];
sx q[0];
rz(1.5623708) q[0];
rz(-pi) q[1];
rz(2.9485547) q[2];
sx q[2];
rz(-1.8503354) q[2];
sx q[2];
rz(-2.2045731) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.78649) q[1];
sx q[1];
rz(-0.10663248) q[1];
sx q[1];
rz(2.998888) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1920131) q[3];
sx q[3];
rz(-2.5887244) q[3];
sx q[3];
rz(-1.0903996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.55591136) q[2];
sx q[2];
rz(-0.53148091) q[2];
sx q[2];
rz(1.4038203) q[2];
rz(-2.3445271) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-1.4306915) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(-2.8523493) q[0];
rz(-0.62492433) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(1.8274868) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02567357) q[0];
sx q[0];
rz(-2.325255) q[0];
sx q[0];
rz(-1.4195819) q[0];
rz(-pi) q[1];
rz(-2.3938789) q[2];
sx q[2];
rz(-1.3249825) q[2];
sx q[2];
rz(-0.73033787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1738759) q[1];
sx q[1];
rz(-1.6248676) q[1];
sx q[1];
rz(1.4500344) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30808361) q[3];
sx q[3];
rz(-2.3039673) q[3];
sx q[3];
rz(-1.777491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73359314) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(-1.7129664) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6190417) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(1.8956986) q[0];
rz(-3.030581) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(2.5949809) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7738757) q[0];
sx q[0];
rz(-1.1989294) q[0];
sx q[0];
rz(-1.0982151) q[0];
x q[1];
rz(2.390929) q[2];
sx q[2];
rz(-2.6235136) q[2];
sx q[2];
rz(-0.93875611) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6706898) q[1];
sx q[1];
rz(-2.7109475) q[1];
sx q[1];
rz(1.9914658) q[1];
rz(-pi) q[2];
rz(-0.9162174) q[3];
sx q[3];
rz(-1.9700053) q[3];
sx q[3];
rz(0.44772128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0299915) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(-1.6938422) q[2];
rz(-1.1374121) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(0.64731961) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2820213) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(2.4243673) q[0];
rz(1.9316797) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-2.4338914) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4004138) q[0];
sx q[0];
rz(-0.78010633) q[0];
sx q[0];
rz(1.0723423) q[0];
rz(-2.999448) q[2];
sx q[2];
rz(-1.3488349) q[2];
sx q[2];
rz(2.1665426) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51703875) q[1];
sx q[1];
rz(-2.3561764) q[1];
sx q[1];
rz(2.024827) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35499771) q[3];
sx q[3];
rz(-2.1600351) q[3];
sx q[3];
rz(2.0139351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6282965) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(0.90325242) q[2];
rz(1.5385657) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(-0.8159591) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(-2.6782398) q[2];
sx q[2];
rz(-1.9909161) q[2];
sx q[2];
rz(-1.9009895) q[2];
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
