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
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(-1.0645359) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2250741) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(-2.201447) q[0];
x q[1];
rz(2.0374374) q[2];
sx q[2];
rz(-2.4541514) q[2];
sx q[2];
rz(1.3165557) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2805466) q[1];
sx q[1];
rz(-1.4006873) q[1];
sx q[1];
rz(2.4982846) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75818054) q[3];
sx q[3];
rz(-1.7054134) q[3];
sx q[3];
rz(1.7077703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(2.1814363) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(-1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(0.91645855) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(0.8786456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71683305) q[0];
sx q[0];
rz(-0.75100198) q[0];
sx q[0];
rz(-0.99320937) q[0];
x q[1];
rz(0.61848817) q[2];
sx q[2];
rz(-1.3534708) q[2];
sx q[2];
rz(2.6478812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0658873) q[1];
sx q[1];
rz(-1.9888708) q[1];
sx q[1];
rz(-1.651152) q[1];
rz(-1.2611748) q[3];
sx q[3];
rz(-1.505758) q[3];
sx q[3];
rz(1.0950973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8615222) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(0.64727616) q[2];
rz(-0.17368008) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(-1.9672111) q[1];
sx q[1];
rz(-2.0203967) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5593528) q[0];
sx q[0];
rz(-1.9358578) q[0];
sx q[0];
rz(2.3854726) q[0];
x q[1];
rz(-1.7240702) q[2];
sx q[2];
rz(-0.3970662) q[2];
sx q[2];
rz(-3.0561471) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.725863) q[1];
sx q[1];
rz(-1.2408857) q[1];
sx q[1];
rz(1.2582474) q[1];
x q[2];
rz(-0.49384533) q[3];
sx q[3];
rz(-1.4673127) q[3];
sx q[3];
rz(-2.8007357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(0.95139727) q[2];
rz(-2.4915063) q[3];
sx q[3];
rz(-1.254046) q[3];
sx q[3];
rz(-2.8459809) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6275416) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(-0.78805584) q[0];
rz(-2.0846941) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-0.11638164) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48556337) q[0];
sx q[0];
rz(-0.93755994) q[0];
sx q[0];
rz(-1.602519) q[0];
rz(-2.1999947) q[2];
sx q[2];
rz(-0.59855748) q[2];
sx q[2];
rz(-1.4231921) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.14028206) q[1];
sx q[1];
rz(-1.9754859) q[1];
sx q[1];
rz(-0.54108041) q[1];
rz(-pi) q[2];
rz(1.6223525) q[3];
sx q[3];
rz(-2.0600852) q[3];
sx q[3];
rz(1.4435022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(2.8895565) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(0.56548059) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(-0.53260032) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.1486357) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2411597) q[0];
sx q[0];
rz(-2.1923089) q[0];
sx q[0];
rz(3.050699) q[0];
rz(-0.45102851) q[2];
sx q[2];
rz(-2.1941059) q[2];
sx q[2];
rz(-0.23652467) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3061805) q[1];
sx q[1];
rz(-2.7369943) q[1];
sx q[1];
rz(2.6447891) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7025181) q[3];
sx q[3];
rz(-0.77507654) q[3];
sx q[3];
rz(-2.4174158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1308412) q[2];
sx q[2];
rz(-2.4804219) q[2];
sx q[2];
rz(2.1095236) q[2];
rz(-2.4268835) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-2.9352303) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26158953) q[0];
sx q[0];
rz(-0.82669175) q[0];
sx q[0];
rz(-1.6893301) q[0];
x q[1];
rz(0.66531078) q[2];
sx q[2];
rz(-2.1449001) q[2];
sx q[2];
rz(-2.2523508) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1604662) q[1];
sx q[1];
rz(-1.4741352) q[1];
sx q[1];
rz(1.4658982) q[1];
rz(1.2431074) q[3];
sx q[3];
rz(-1.7565691) q[3];
sx q[3];
rz(-2.1873308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53987327) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(-1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(-1.9304088) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9528708) q[0];
sx q[0];
rz(-2.0824008) q[0];
sx q[0];
rz(-2.1923724) q[0];
rz(-pi) q[1];
rz(-0.25522916) q[2];
sx q[2];
rz(-2.6380939) q[2];
sx q[2];
rz(-0.90571678) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32528977) q[1];
sx q[1];
rz(-1.6386697) q[1];
sx q[1];
rz(1.7273278) q[1];
x q[2];
rz(0.96774775) q[3];
sx q[3];
rz(-1.839404) q[3];
sx q[3];
rz(0.3097765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3367735) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(3.1170735) q[2];
rz(2.426614) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(2.1210282) q[0];
rz(2.9868946) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.221009) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3295791) q[0];
sx q[0];
rz(-2.3276969) q[0];
sx q[0];
rz(-1.9474067) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9646655) q[2];
sx q[2];
rz(-1.7311586) q[2];
sx q[2];
rz(2.4307876) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.11215969) q[1];
sx q[1];
rz(-0.1754079) q[1];
sx q[1];
rz(0.68316858) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8052748) q[3];
sx q[3];
rz(-0.85018966) q[3];
sx q[3];
rz(-1.1433126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(0.60338902) q[3];
sx q[3];
rz(-2.2796977) q[3];
sx q[3];
rz(1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0309546) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(3.0045793) q[0];
rz(2.5367472) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(-0.12577122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14995689) q[0];
sx q[0];
rz(-1.7466674) q[0];
sx q[0];
rz(-1.3791023) q[0];
rz(-pi) q[1];
rz(1.3734829) q[2];
sx q[2];
rz(-1.8981877) q[2];
sx q[2];
rz(1.9257853) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3198493) q[1];
sx q[1];
rz(-2.2935767) q[1];
sx q[1];
rz(0.38139947) q[1];
rz(-pi) q[2];
rz(-0.29571663) q[3];
sx q[3];
rz(-1.3774646) q[3];
sx q[3];
rz(-2.8545692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13828364) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(-0.62018958) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(1.2003157) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(-2.0013924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32894293) q[0];
sx q[0];
rz(-2.86033) q[0];
sx q[0];
rz(-1.2055231) q[0];
rz(-pi) q[1];
x q[1];
rz(2.033038) q[2];
sx q[2];
rz(-1.3369563) q[2];
sx q[2];
rz(0.086694593) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4789341) q[1];
sx q[1];
rz(-1.7591898) q[1];
sx q[1];
rz(1.4807832) q[1];
rz(1.8885918) q[3];
sx q[3];
rz(-0.80298775) q[3];
sx q[3];
rz(3.0190937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7662979) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(2.2429788) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469149) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(0.5207516) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(0.49171369) q[2];
sx q[2];
rz(-1.9946873) q[2];
sx q[2];
rz(-0.36110525) q[2];
rz(0.58581523) q[3];
sx q[3];
rz(-1.2231493) q[3];
sx q[3];
rz(2.3412658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];