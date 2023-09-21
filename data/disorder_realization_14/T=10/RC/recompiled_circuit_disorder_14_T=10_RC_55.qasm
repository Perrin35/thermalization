OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.68552652) q[0];
sx q[0];
rz(-2.752562) q[0];
sx q[0];
rz(-2.2580137) q[0];
rz(3.1318624) q[1];
sx q[1];
rz(4.598773) q[1];
sx q[1];
rz(7.4814319) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0449013) q[0];
sx q[0];
rz(-2.9400819) q[0];
sx q[0];
rz(-2.6411396) q[0];
rz(-0.77550195) q[2];
sx q[2];
rz(-1.8047793) q[2];
sx q[2];
rz(-2.6393059) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8625496) q[1];
sx q[1];
rz(-2.0837796) q[1];
sx q[1];
rz(1.951745) q[1];
x q[2];
rz(-0.046273307) q[3];
sx q[3];
rz(-1.4001906) q[3];
sx q[3];
rz(1.5266649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9464232) q[2];
sx q[2];
rz(-2.158458) q[2];
sx q[2];
rz(-2.9602489) q[2];
rz(-2.8803853) q[3];
sx q[3];
rz(-1.8758592) q[3];
sx q[3];
rz(-0.75631022) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392035) q[0];
sx q[0];
rz(-0.2897245) q[0];
sx q[0];
rz(2.7547577) q[0];
rz(-2.6392) q[1];
sx q[1];
rz(-2.1680809) q[1];
sx q[1];
rz(-1.5997255) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5889266) q[0];
sx q[0];
rz(-0.98836556) q[0];
sx q[0];
rz(-1.0731339) q[0];
rz(1.2453169) q[2];
sx q[2];
rz(-2.6633334) q[2];
sx q[2];
rz(1.1748199) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.220037) q[1];
sx q[1];
rz(-1.5350635) q[1];
sx q[1];
rz(1.5608556) q[1];
x q[2];
rz(-0.76962556) q[3];
sx q[3];
rz(-2.5105021) q[3];
sx q[3];
rz(-2.1345994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5388422) q[2];
sx q[2];
rz(-2.2637612) q[2];
sx q[2];
rz(-1.4146457) q[2];
rz(2.291262) q[3];
sx q[3];
rz(-0.43262216) q[3];
sx q[3];
rz(1.458228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(0.26281115) q[0];
sx q[0];
rz(-1.5314064) q[0];
sx q[0];
rz(2.5313654) q[0];
rz(1.8521076) q[1];
sx q[1];
rz(-0.97924966) q[1];
sx q[1];
rz(-2.1496225) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45527601) q[0];
sx q[0];
rz(-1.0193829) q[0];
sx q[0];
rz(-2.9247012) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8875071) q[2];
sx q[2];
rz(-0.68426364) q[2];
sx q[2];
rz(-0.64988818) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85769535) q[1];
sx q[1];
rz(-1.1313492) q[1];
sx q[1];
rz(-1.6698014) q[1];
rz(1.3258341) q[3];
sx q[3];
rz(-2.0162597) q[3];
sx q[3];
rz(-0.20416343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.758574) q[2];
sx q[2];
rz(-2.980361) q[2];
sx q[2];
rz(0.26829159) q[2];
rz(0.39408436) q[3];
sx q[3];
rz(-1.2309309) q[3];
sx q[3];
rz(-0.18850732) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85369337) q[0];
sx q[0];
rz(-0.51369602) q[0];
sx q[0];
rz(-2.7365141) q[0];
rz(-2.4515117) q[1];
sx q[1];
rz(-1.1578553) q[1];
sx q[1];
rz(2.4437723) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9245802) q[0];
sx q[0];
rz(-0.092886535) q[0];
sx q[0];
rz(-1.3148944) q[0];
rz(-1.1890829) q[2];
sx q[2];
rz(-1.0587947) q[2];
sx q[2];
rz(2.2005759) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.78742541) q[1];
sx q[1];
rz(-1.0855055) q[1];
sx q[1];
rz(0.20809681) q[1];
rz(-1.1504193) q[3];
sx q[3];
rz(-1.8432518) q[3];
sx q[3];
rz(0.099893173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4219389) q[2];
sx q[2];
rz(-1.3976588) q[2];
sx q[2];
rz(1.6323803) q[2];
rz(0.40431067) q[3];
sx q[3];
rz(-0.68250889) q[3];
sx q[3];
rz(-1.6633165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835835) q[0];
sx q[0];
rz(-1.785935) q[0];
sx q[0];
rz(-2.1160545) q[0];
rz(-0.57199663) q[1];
sx q[1];
rz(-1.0943741) q[1];
sx q[1];
rz(-0.62932032) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9176661) q[0];
sx q[0];
rz(-1.9804945) q[0];
sx q[0];
rz(1.8976589) q[0];
x q[1];
rz(-1.8848558) q[2];
sx q[2];
rz(-0.70054189) q[2];
sx q[2];
rz(2.8741921) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8992936) q[1];
sx q[1];
rz(-1.0596024) q[1];
sx q[1];
rz(-1.7023229) q[1];
rz(-2.1719645) q[3];
sx q[3];
rz(-1.3519577) q[3];
sx q[3];
rz(2.9588483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5489674) q[2];
sx q[2];
rz(-2.7719345) q[2];
sx q[2];
rz(-0.34234753) q[2];
rz(1.7957548) q[3];
sx q[3];
rz(-1.4474409) q[3];
sx q[3];
rz(2.9455744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65258566) q[0];
sx q[0];
rz(-1.2359897) q[0];
sx q[0];
rz(-2.956399) q[0];
rz(1.406503) q[1];
sx q[1];
rz(-1.0909189) q[1];
sx q[1];
rz(-1.3669744) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0104116) q[0];
sx q[0];
rz(-2.512012) q[0];
sx q[0];
rz(-2.7601943) q[0];
x q[1];
rz(-2.30079) q[2];
sx q[2];
rz(-1.3422988) q[2];
sx q[2];
rz(-2.218354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5806611) q[1];
sx q[1];
rz(-1.4189737) q[1];
sx q[1];
rz(3.0488357) q[1];
x q[2];
rz(-1.1432511) q[3];
sx q[3];
rz(-1.0435836) q[3];
sx q[3];
rz(0.15641071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.24484816) q[2];
sx q[2];
rz(-2.5881793) q[2];
sx q[2];
rz(2.2582167) q[2];
rz(-1.4128489) q[3];
sx q[3];
rz(-2.4491375) q[3];
sx q[3];
rz(-2.3278055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8900523) q[0];
sx q[0];
rz(-0.10245704) q[0];
sx q[0];
rz(1.2782156) q[0];
rz(-0.037840769) q[1];
sx q[1];
rz(-2.3262639) q[1];
sx q[1];
rz(1.3758804) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28988518) q[0];
sx q[0];
rz(-0.87806784) q[0];
sx q[0];
rz(-0.38498199) q[0];
rz(-pi) q[1];
rz(0.059869754) q[2];
sx q[2];
rz(-1.0979568) q[2];
sx q[2];
rz(1.5945895) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3657276) q[1];
sx q[1];
rz(-0.15656808) q[1];
sx q[1];
rz(-2.4003029) q[1];
rz(-pi) q[2];
rz(1.901058) q[3];
sx q[3];
rz(-1.5152647) q[3];
sx q[3];
rz(-0.33952573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6841131) q[2];
sx q[2];
rz(-0.71762466) q[2];
sx q[2];
rz(0.73105556) q[2];
rz(3.030792) q[3];
sx q[3];
rz(-1.5564857) q[3];
sx q[3];
rz(-2.4462162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59657997) q[0];
sx q[0];
rz(-0.87485635) q[0];
sx q[0];
rz(2.4080283) q[0];
rz(-2.5336174) q[1];
sx q[1];
rz(-1.9476451) q[1];
sx q[1];
rz(0.2342934) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83345862) q[0];
sx q[0];
rz(-1.8292384) q[0];
sx q[0];
rz(1.0779557) q[0];
x q[1];
rz(0.73306002) q[2];
sx q[2];
rz(-0.86793938) q[2];
sx q[2];
rz(1.5066063) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1842321) q[1];
sx q[1];
rz(-2.3174274) q[1];
sx q[1];
rz(-2.8824473) q[1];
rz(-pi) q[2];
rz(-0.33741823) q[3];
sx q[3];
rz(-2.1302345) q[3];
sx q[3];
rz(0.19677256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0155448) q[2];
sx q[2];
rz(-1.8323106) q[2];
sx q[2];
rz(1.4253433) q[2];
rz(1.4632633) q[3];
sx q[3];
rz(-2.3571456) q[3];
sx q[3];
rz(1.6459758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54206806) q[0];
sx q[0];
rz(-0.33518377) q[0];
sx q[0];
rz(-1.9482127) q[0];
rz(1.880973) q[1];
sx q[1];
rz(-1.7766989) q[1];
sx q[1];
rz(-0.9448005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.955907) q[0];
sx q[0];
rz(-0.80230306) q[0];
sx q[0];
rz(2.5923652) q[0];
rz(-pi) q[1];
rz(2.2009732) q[2];
sx q[2];
rz(-0.47626469) q[2];
sx q[2];
rz(2.1911052) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33465696) q[1];
sx q[1];
rz(-1.2922704) q[1];
sx q[1];
rz(1.7486228) q[1];
rz(-2.0637022) q[3];
sx q[3];
rz(-0.3993881) q[3];
sx q[3];
rz(-0.82769094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21645674) q[2];
sx q[2];
rz(-2.0704806) q[2];
sx q[2];
rz(-0.28820583) q[2];
rz(0.47973412) q[3];
sx q[3];
rz(-1.0498485) q[3];
sx q[3];
rz(1.5079927) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11186803) q[0];
sx q[0];
rz(-0.27619633) q[0];
sx q[0];
rz(0.9129886) q[0];
rz(0.37462014) q[1];
sx q[1];
rz(-1.4034142) q[1];
sx q[1];
rz(0.8909117) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83715314) q[0];
sx q[0];
rz(-1.8189948) q[0];
sx q[0];
rz(-2.7960143) q[0];
rz(2.183379) q[2];
sx q[2];
rz(-1.805086) q[2];
sx q[2];
rz(-0.9226375) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3112463) q[1];
sx q[1];
rz(-1.3724569) q[1];
sx q[1];
rz(2.7970008) q[1];
x q[2];
rz(2.9534146) q[3];
sx q[3];
rz(-1.6450226) q[3];
sx q[3];
rz(1.0159514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9051819) q[2];
sx q[2];
rz(-0.85835251) q[2];
sx q[2];
rz(-0.50160828) q[2];
rz(1.8524528) q[3];
sx q[3];
rz(-1.4533549) q[3];
sx q[3];
rz(2.6954209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63507737) q[0];
sx q[0];
rz(-1.4415393) q[0];
sx q[0];
rz(-2.517979) q[0];
rz(-2.0093909) q[1];
sx q[1];
rz(-0.75695801) q[1];
sx q[1];
rz(-3.0523041) q[1];
rz(-0.51433993) q[2];
sx q[2];
rz(-0.30403501) q[2];
sx q[2];
rz(-2.4652849) q[2];
rz(-2.1174259) q[3];
sx q[3];
rz(-1.5285986) q[3];
sx q[3];
rz(-0.91615208) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
