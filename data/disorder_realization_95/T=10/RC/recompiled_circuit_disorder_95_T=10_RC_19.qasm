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
rz(2.907213) q[1];
sx q[1];
rz(-2.8657764) q[1];
sx q[1];
rz(-2.0770567) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9165186) q[0];
sx q[0];
rz(-1.6406035) q[0];
sx q[0];
rz(-0.9401456) q[0];
x q[1];
rz(1.1041553) q[2];
sx q[2];
rz(-0.68744126) q[2];
sx q[2];
rz(-1.825037) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.068163888) q[1];
sx q[1];
rz(-2.4792719) q[1];
sx q[1];
rz(2.8627002) q[1];
x q[2];
rz(0.1944794) q[3];
sx q[3];
rz(-2.3739061) q[3];
sx q[3];
rz(3.1377813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(-0.73195362) q[2];
rz(-0.96015635) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17523781) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(-2.6610999) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(0.8786456) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40936138) q[0];
sx q[0];
rz(-1.9525813) q[0];
sx q[0];
rz(2.2344927) q[0];
rz(-pi) q[1];
rz(-2.5231045) q[2];
sx q[2];
rz(-1.7881219) q[2];
sx q[2];
rz(-2.6478812) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2616927) q[1];
sx q[1];
rz(-0.42527929) q[1];
sx q[1];
rz(2.9628423) q[1];
rz(-pi) q[2];
rz(1.7813803) q[3];
sx q[3];
rz(-0.31616351) q[3];
sx q[3];
rz(2.8663243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(0.64727616) q[2];
rz(2.9679126) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(2.98996) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52755255) q[0];
sx q[0];
rz(-2.5019167) q[0];
sx q[0];
rz(-0.24599427) q[0];
rz(-1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(-1.1211959) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778053) q[0];
sx q[0];
rz(-2.2664547) q[0];
sx q[0];
rz(-2.0545161) q[0];
rz(-pi) q[1];
rz(1.9636743) q[2];
sx q[2];
rz(-1.5117206) q[2];
sx q[2];
rz(-1.7977561) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0507625) q[1];
sx q[1];
rz(-1.275626) q[1];
sx q[1];
rz(2.7961618) q[1];
rz(-2.6477473) q[3];
sx q[3];
rz(-1.6742799) q[3];
sx q[3];
rz(-2.8007357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-2.1901954) q[2];
rz(0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
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
rz(3.025211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1040092) q[0];
sx q[0];
rz(-1.5452256) q[0];
sx q[0];
rz(0.63347647) q[0];
rz(-pi) q[1];
rz(2.0747244) q[2];
sx q[2];
rz(-1.2328086) q[2];
sx q[2];
rz(-0.68901125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.14028206) q[1];
sx q[1];
rz(-1.1661068) q[1];
sx q[1];
rz(-0.54108041) q[1];
rz(1.6223525) q[3];
sx q[3];
rz(-2.0600852) q[3];
sx q[3];
rz(1.4435022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(-0.25203618) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-2.9791322) q[3];
sx q[3];
rz(-2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(1.4415007) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(1.1486357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7589353) q[0];
sx q[0];
rz(-1.6446582) q[0];
sx q[0];
rz(0.9473243) q[0];
x q[1];
rz(1.0257452) q[2];
sx q[2];
rz(-2.3902241) q[2];
sx q[2];
rz(0.92912208) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.83541218) q[1];
sx q[1];
rz(-0.40459834) q[1];
sx q[1];
rz(-2.6447891) q[1];
rz(-pi) q[2];
rz(0.72539056) q[3];
sx q[3];
rz(-1.2687506) q[3];
sx q[3];
rz(-2.6186752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0107515) q[2];
sx q[2];
rz(-2.4804219) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(2.7456039) q[0];
rz(-1.4453325) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(0.2063624) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9128742) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(2.3939783) q[0];
rz(-pi) q[1];
rz(0.88271898) q[2];
sx q[2];
rz(-2.1157017) q[2];
sx q[2];
rz(2.8628652) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3316162) q[1];
sx q[1];
rz(-2.9990701) q[1];
sx q[1];
rz(-2.317821) q[1];
rz(-pi) q[2];
rz(1.0422802) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(-0.11881766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53987327) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(-0.27077857) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(1.9304088) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1887218) q[0];
sx q[0];
rz(-1.0591918) q[0];
sx q[0];
rz(2.1923724) q[0];
rz(0.48970512) q[2];
sx q[2];
rz(-1.4486794) q[2];
sx q[2];
rz(0.44039886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8163029) q[1];
sx q[1];
rz(-1.5029229) q[1];
sx q[1];
rz(-1.4142649) q[1];
rz(-2.8190523) q[3];
sx q[3];
rz(-2.1493559) q[3];
sx q[3];
rz(-1.4417779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(-0.024519196) q[2];
rz(-0.71497861) q[3];
sx q[3];
rz(-1.4638126) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0054758469) q[0];
sx q[0];
rz(-1.5908717) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(-2.9868946) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.9205836) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1655501) q[0];
sx q[0];
rz(-1.8414458) q[0];
sx q[0];
rz(0.793215) q[0];
rz(-pi) q[1];
rz(2.3983725) q[2];
sx q[2];
rz(-0.23822242) q[2];
sx q[2];
rz(-1.5889578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80291058) q[1];
sx q[1];
rz(-1.435034) q[1];
sx q[1];
rz(1.6822097) q[1];
x q[2];
rz(-0.33631781) q[3];
sx q[3];
rz(-0.85018966) q[3];
sx q[3];
rz(1.99828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(2.4364046) q[2];
rz(0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0309546) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(-0.13701339) q[0];
rz(-0.6048454) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(-0.12577122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1544979) q[0];
sx q[0];
rz(-2.8821766) q[0];
sx q[0];
rz(2.3214066) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33340402) q[2];
sx q[2];
rz(-1.3840904) q[2];
sx q[2];
rz(-2.8508027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3198493) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(-2.7601932) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59154193) q[3];
sx q[3];
rz(-0.35174832) q[3];
sx q[3];
rz(-1.2951617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13828364) q[2];
sx q[2];
rz(-2.1676899) q[2];
sx q[2];
rz(-0.70927817) q[2];
rz(0.62018958) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1518635) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(1.9412769) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(2.0013924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88975554) q[0];
sx q[0];
rz(-1.4714843) q[0];
sx q[0];
rz(-1.8343783) q[0];
rz(-1.0802202) q[2];
sx q[2];
rz(-0.51418257) q[2];
sx q[2];
rz(-1.0487923) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1117489) q[1];
sx q[1];
rz(-2.9330301) q[1];
sx q[1];
rz(-2.7010121) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31302281) q[3];
sx q[3];
rz(-0.81837294) q[3];
sx q[3];
rz(-0.5648053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-2.2429788) q[2];
rz(-3.1344154) q[3];
sx q[3];
rz(-2.3982748) q[3];
sx q[3];
rz(-1.3557419) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2469149) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.649879) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(-1.1605916) q[3];
sx q[3];
rz(-2.1174178) q[3];
sx q[3];
rz(-2.5934364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
