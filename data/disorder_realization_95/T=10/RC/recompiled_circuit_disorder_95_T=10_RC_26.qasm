OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(-1.002123) q[0];
sx q[0];
rz(2.2440417) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(2.0770567) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8911154) q[0];
sx q[0];
rz(-0.63397898) q[0];
sx q[0];
rz(1.4527713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93809442) q[2];
sx q[2];
rz(-1.860306) q[2];
sx q[2];
rz(3.0245568) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2805466) q[1];
sx q[1];
rz(-1.7409054) q[1];
sx q[1];
rz(2.4982846) q[1];
rz(-0.75818054) q[3];
sx q[3];
rz(-1.7054134) q[3];
sx q[3];
rz(1.4338223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.66427461) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(2.409639) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(0.91645855) q[0];
rz(-0.48049277) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(-2.2629471) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4449578) q[0];
sx q[0];
rz(-2.1793471) q[0];
sx q[0];
rz(0.47136013) q[0];
rz(2.5231045) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(-2.6478812) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.87989992) q[1];
sx q[1];
rz(-2.7163134) q[1];
sx q[1];
rz(-2.9628423) q[1];
rz(-pi) q[2];
rz(-1.2611748) q[3];
sx q[3];
rz(-1.6358346) q[3];
sx q[3];
rz(-1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-0.64727616) q[2];
rz(0.17368008) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52755255) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(-1.4100769) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(2.0203967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5593528) q[0];
sx q[0];
rz(-1.9358578) q[0];
sx q[0];
rz(-2.3854726) q[0];
rz(0.063935117) q[2];
sx q[2];
rz(-1.1786412) q[2];
sx q[2];
rz(0.25142297) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9417291) q[1];
sx q[1];
rz(-0.45048303) q[1];
sx q[1];
rz(0.73168879) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9259053) q[3];
sx q[3];
rz(-0.50369278) q[3];
sx q[3];
rz(1.7220725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-2.1901954) q[2];
rz(-2.4915063) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(-2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(0.78805584) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(0.11638164) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48556337) q[0];
sx q[0];
rz(-0.93755994) q[0];
sx q[0];
rz(-1.602519) q[0];
rz(0.38168455) q[2];
sx q[2];
rz(-2.0437864) q[2];
sx q[2];
rz(0.7009398) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1981922) q[1];
sx q[1];
rz(-1.0775837) q[1];
sx q[1];
rz(-2.0342159) q[1];
rz(-pi) q[2];
rz(1.5192401) q[3];
sx q[3];
rz(-2.0600852) q[3];
sx q[3];
rz(1.6980905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(-0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(-2.6089923) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.9929569) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2411597) q[0];
sx q[0];
rz(-0.94928375) q[0];
sx q[0];
rz(0.09089367) q[0];
x q[1];
rz(-1.0257452) q[2];
sx q[2];
rz(-2.3902241) q[2];
sx q[2];
rz(2.2124706) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3061805) q[1];
sx q[1];
rz(-0.40459834) q[1];
sx q[1];
rz(0.49680357) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4162021) q[3];
sx q[3];
rz(-1.8728421) q[3];
sx q[3];
rz(0.52291742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0107515) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(-0.71470913) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5500568) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(-0.39598879) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(0.2063624) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26158953) q[0];
sx q[0];
rz(-0.82669175) q[0];
sx q[0];
rz(1.6893301) q[0];
rz(-0.66531078) q[2];
sx q[2];
rz(-2.1449001) q[2];
sx q[2];
rz(-0.88924185) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.57950912) q[1];
sx q[1];
rz(-1.4663896) q[1];
sx q[1];
rz(0.097192055) q[1];
rz(-pi) q[2];
rz(2.9456375) q[3];
sx q[3];
rz(-1.8926419) q[3];
sx q[3];
rz(-2.4623507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6017194) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(2.3932636) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(1.07871) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(2.5653429) q[0];
rz(2.9684864) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(-1.2111838) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9241087) q[0];
sx q[0];
rz(-2.3586914) q[0];
sx q[0];
rz(2.3379675) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48970512) q[2];
sx q[2];
rz(-1.6929132) q[2];
sx q[2];
rz(0.44039886) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8853828) q[1];
sx q[1];
rz(-1.7269644) q[1];
sx q[1];
rz(-0.068710879) q[1];
rz(-0.96774775) q[3];
sx q[3];
rz(-1.3021886) q[3];
sx q[3];
rz(-2.8318162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(0.024519196) q[2];
rz(-0.71497861) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(-1.0205644) q[0];
rz(-2.9868946) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(-1.9205836) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8120136) q[0];
sx q[0];
rz(-2.3276969) q[0];
sx q[0];
rz(-1.9474067) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74322015) q[2];
sx q[2];
rz(-2.9033702) q[2];
sx q[2];
rz(1.5889578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.358566) q[1];
sx q[1];
rz(-1.6811803) q[1];
sx q[1];
rz(-0.13659887) q[1];
rz(2.8052748) q[3];
sx q[3];
rz(-0.85018966) q[3];
sx q[3];
rz(1.99828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.28875479) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(2.4364046) q[2];
rz(0.60338902) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(-1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.110638) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(-3.0045793) q[0];
rz(-2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(-3.0158214) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14995689) q[0];
sx q[0];
rz(-1.7466674) q[0];
sx q[0];
rz(-1.3791023) q[0];
x q[1];
rz(2.618082) q[2];
sx q[2];
rz(-0.38041174) q[2];
sx q[2];
rz(1.3695804) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3198493) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(2.7601932) q[1];
rz(-pi) q[2];
x q[2];
rz(2.845876) q[3];
sx q[3];
rz(-1.7641281) q[3];
sx q[3];
rz(-0.28702345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13828364) q[2];
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
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(-1.1402003) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2518371) q[0];
sx q[0];
rz(-1.4714843) q[0];
sx q[0];
rz(-1.3072144) q[0];
x q[1];
rz(-1.0802202) q[2];
sx q[2];
rz(-2.6274101) q[2];
sx q[2];
rz(-2.0928004) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1117489) q[1];
sx q[1];
rz(-0.20856253) q[1];
sx q[1];
rz(2.7010121) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8885918) q[3];
sx q[3];
rz(-0.80298775) q[3];
sx q[3];
rz(-3.0190937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(0.89861384) q[2];
rz(-0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89467775) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(-2.3788135) q[2];
sx q[2];
rz(-2.503958) q[2];
sx q[2];
rz(-2.5867953) q[2];
rz(1.981001) q[3];
sx q[3];
rz(-2.1174178) q[3];
sx q[3];
rz(-2.5934364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
