OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(-2.2440417) q[0];
rz(-3.3759723) q[1];
sx q[1];
rz(3.4174089) q[1];
sx q[1];
rz(13.630907) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8911154) q[0];
sx q[0];
rz(-0.63397898) q[0];
sx q[0];
rz(-1.6888213) q[0];
x q[1];
rz(0.35383309) q[2];
sx q[2];
rz(-2.1733123) q[2];
sx q[2];
rz(1.2474071) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2805466) q[1];
sx q[1];
rz(-1.4006873) q[1];
sx q[1];
rz(-0.64330805) q[1];
rz(1.7552056) q[3];
sx q[3];
rz(-2.3204436) q[3];
sx q[3];
rz(2.8781995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66427461) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.9663548) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(0.8786456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7322313) q[0];
sx q[0];
rz(-1.9525813) q[0];
sx q[0];
rz(-0.9071) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61848817) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(-2.6478812) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2616927) q[1];
sx q[1];
rz(-0.42527929) q[1];
sx q[1];
rz(-2.9628423) q[1];
rz(1.2611748) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(-0.64727616) q[2];
rz(-0.17368008) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(-2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6140401) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(-1.4100769) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-1.1211959) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5593528) q[0];
sx q[0];
rz(-1.2057349) q[0];
sx q[0];
rz(2.3854726) q[0];
rz(-1.7240702) q[2];
sx q[2];
rz(-0.3970662) q[2];
sx q[2];
rz(0.085445554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0908302) q[1];
sx q[1];
rz(-1.275626) q[1];
sx q[1];
rz(-2.7961618) q[1];
rz(0.21568732) q[3];
sx q[3];
rz(-2.6378999) q[3];
sx q[3];
rz(-1.4195201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7401509) q[2];
sx q[2];
rz(-1.695305) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(2.3535368) q[0];
rz(2.0846941) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(3.025211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43198904) q[0];
sx q[0];
rz(-2.5076712) q[0];
sx q[0];
rz(3.0984127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.94159796) q[2];
sx q[2];
rz(-2.5430352) q[2];
sx q[2];
rz(1.7184005) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9434005) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(1.1073768) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0451123) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(1.5887367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(-0.37825545) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-2.7799515) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(2.6089923) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(-1.1486357) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.085885) q[0];
sx q[0];
rz(-0.62725337) q[0];
sx q[0];
rz(-1.4447312) q[0];
rz(-2.1158475) q[2];
sx q[2];
rz(-2.3902241) q[2];
sx q[2];
rz(-2.2124706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3682813) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(1.3694622) q[1];
x q[2];
rz(1.1762189) q[3];
sx q[3];
rz(-0.88486457) q[3];
sx q[3];
rz(-1.3057614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1308412) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(-2.1095236) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-0.023795279) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5500568) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(2.7456039) q[0];
rz(-1.4453325) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-2.9352303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43564046) q[0];
sx q[0];
rz(-2.3898976) q[0];
sx q[0];
rz(-3.0138426) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2588737) q[2];
sx q[2];
rz(-1.025891) q[2];
sx q[2];
rz(-2.8628652) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.57950912) q[1];
sx q[1];
rz(-1.675203) q[1];
sx q[1];
rz(-3.0444006) q[1];
rz(-2.0993125) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(-0.11881766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53987327) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(-0.27077857) q[2];
rz(-2.3932636) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2449743) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(-2.9684864) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(-1.2111838) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9241087) q[0];
sx q[0];
rz(-2.3586914) q[0];
sx q[0];
rz(0.80362513) q[0];
x q[1];
rz(-0.48970512) q[2];
sx q[2];
rz(-1.4486794) q[2];
sx q[2];
rz(2.7011938) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8853828) q[1];
sx q[1];
rz(-1.7269644) q[1];
sx q[1];
rz(-3.0728818) q[1];
rz(-pi) q[2];
rz(2.022642) q[3];
sx q[3];
rz(-2.4882836) q[3];
sx q[3];
rz(-2.2484231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(-0.024519196) q[2];
rz(-2.426614) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(1.0205644) q[0];
rz(-2.9868946) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(-1.9205836) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1655501) q[0];
sx q[0];
rz(-1.3001469) q[0];
sx q[0];
rz(-2.3483777) q[0];
rz(-pi) q[1];
rz(0.17692716) q[2];
sx q[2];
rz(-1.7311586) q[2];
sx q[2];
rz(-2.4307876) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11215969) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(0.68316858) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3200795) q[3];
sx q[3];
rz(-1.320208) q[3];
sx q[3];
rz(0.65419765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(-2.4364046) q[2];
rz(-2.5382036) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(-1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(-3.0045793) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(0.12577122) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3868956) q[0];
sx q[0];
rz(-1.3820952) q[0];
sx q[0];
rz(0.17908355) q[0];
rz(2.8081886) q[2];
sx q[2];
rz(-1.7575022) q[2];
sx q[2];
rz(-0.29078996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77482624) q[1];
sx q[1];
rz(-0.80087304) q[1];
sx q[1];
rz(-1.1714539) q[1];
rz(1.7726692) q[3];
sx q[3];
rz(-1.8608421) q[3];
sx q[3];
rz(1.2253075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.003309) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(0.70927817) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-1.1518635) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(-2.0013924) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8126497) q[0];
sx q[0];
rz(-2.86033) q[0];
sx q[0];
rz(1.2055231) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8814949) q[2];
sx q[2];
rz(-2.0195228) q[2];
sx q[2];
rz(-1.5425494) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0298437) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(-0.44058056) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79348989) q[3];
sx q[3];
rz(-1.3440545) q[3];
sx q[3];
rz(1.2236809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7662979) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89467775) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-1.0976787) q[2];
sx q[2];
rz(-1.1259148) q[2];
sx q[2];
rz(-1.7150707) q[2];
rz(-1.981001) q[3];
sx q[3];
rz(-1.0241749) q[3];
sx q[3];
rz(0.54815626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];