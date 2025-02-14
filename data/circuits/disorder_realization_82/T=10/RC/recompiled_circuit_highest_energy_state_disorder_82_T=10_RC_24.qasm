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
rz(-3.1373625) q[0];
sx q[0];
rz(-2.5047996) q[0];
sx q[0];
rz(1.1871583) q[0];
rz(4.123426) q[1];
sx q[1];
rz(0.47458664) q[1];
sx q[1];
rz(8.4835806) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6146916) q[0];
sx q[0];
rz(-2.6462835) q[0];
sx q[0];
rz(0.0035347819) q[0];
rz(-pi) q[1];
rz(-1.0082939) q[2];
sx q[2];
rz(-0.23782158) q[2];
sx q[2];
rz(0.49006762) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8601712) q[1];
sx q[1];
rz(-2.1263391) q[1];
sx q[1];
rz(-1.9940328) q[1];
rz(2.0273846) q[3];
sx q[3];
rz(-1.5043882) q[3];
sx q[3];
rz(-0.6765863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5473951) q[2];
sx q[2];
rz(-1.6387458) q[2];
sx q[2];
rz(-1.0198062) q[2];
rz(-2.8504573) q[3];
sx q[3];
rz(-2.9499493) q[3];
sx q[3];
rz(-1.3297133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52784598) q[0];
sx q[0];
rz(-0.83714038) q[0];
sx q[0];
rz(1.5035195) q[0];
rz(-1.9095406) q[1];
sx q[1];
rz(-2.3161395) q[1];
sx q[1];
rz(1.8881316) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37663662) q[0];
sx q[0];
rz(-0.47237637) q[0];
sx q[0];
rz(1.0563904) q[0];
x q[1];
rz(-0.49053662) q[2];
sx q[2];
rz(-0.56714688) q[2];
sx q[2];
rz(0.74037742) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.277569) q[1];
sx q[1];
rz(-1.4269842) q[1];
sx q[1];
rz(-0.20419232) q[1];
rz(-2.0157865) q[3];
sx q[3];
rz(-2.2644694) q[3];
sx q[3];
rz(1.3563658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2848009) q[2];
sx q[2];
rz(-2.5992726) q[2];
sx q[2];
rz(0.33033237) q[2];
rz(-0.10279837) q[3];
sx q[3];
rz(-1.8967352) q[3];
sx q[3];
rz(0.61906329) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61515808) q[0];
sx q[0];
rz(-1.1410843) q[0];
sx q[0];
rz(-2.0779628) q[0];
rz(3.0377153) q[1];
sx q[1];
rz(-2.3492298) q[1];
sx q[1];
rz(2.0361384) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.605382) q[0];
sx q[0];
rz(-0.057312556) q[0];
sx q[0];
rz(-0.97672279) q[0];
rz(-1.8806522) q[2];
sx q[2];
rz(-0.83854691) q[2];
sx q[2];
rz(-2.1149879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6686386) q[1];
sx q[1];
rz(-2.611043) q[1];
sx q[1];
rz(1.3282685) q[1];
rz(0.01171578) q[3];
sx q[3];
rz(-1.612631) q[3];
sx q[3];
rz(1.5277629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1198472) q[2];
sx q[2];
rz(-0.80588078) q[2];
sx q[2];
rz(2.9160807) q[2];
rz(-1.3269904) q[3];
sx q[3];
rz(-0.83695379) q[3];
sx q[3];
rz(-0.64083159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2284112) q[0];
sx q[0];
rz(-1.5629733) q[0];
sx q[0];
rz(-2.5033503) q[0];
rz(-1.111521) q[1];
sx q[1];
rz(-0.94739729) q[1];
sx q[1];
rz(-0.18724719) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.603724) q[0];
sx q[0];
rz(-0.6938254) q[0];
sx q[0];
rz(-2.0286454) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6333605) q[2];
sx q[2];
rz(-0.79088941) q[2];
sx q[2];
rz(-3.084201) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.9065215) q[1];
sx q[1];
rz(-1.491761) q[1];
sx q[1];
rz(2.6071878) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4832014) q[3];
sx q[3];
rz(-0.85932743) q[3];
sx q[3];
rz(-0.30184823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.036666544) q[2];
sx q[2];
rz(-1.5035524) q[2];
sx q[2];
rz(-0.46910134) q[2];
rz(2.8406298) q[3];
sx q[3];
rz(-0.260869) q[3];
sx q[3];
rz(-1.652396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8982573) q[0];
sx q[0];
rz(-1.6247592) q[0];
sx q[0];
rz(0.54018706) q[0];
rz(-0.46145269) q[1];
sx q[1];
rz(-1.6199473) q[1];
sx q[1];
rz(0.25925055) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.777939) q[0];
sx q[0];
rz(-2.724132) q[0];
sx q[0];
rz(2.4221942) q[0];
rz(2.2770693) q[2];
sx q[2];
rz(-0.63535035) q[2];
sx q[2];
rz(0.29624789) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.66189761) q[1];
sx q[1];
rz(-2.5475218) q[1];
sx q[1];
rz(1.8516783) q[1];
x q[2];
rz(-2.2336273) q[3];
sx q[3];
rz(-2.7304711) q[3];
sx q[3];
rz(-3.079252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.79869142) q[2];
sx q[2];
rz(-2.4553802) q[2];
sx q[2];
rz(-0.9101103) q[2];
rz(2.5765007) q[3];
sx q[3];
rz(-1.3684973) q[3];
sx q[3];
rz(-2.4250987) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703576) q[0];
sx q[0];
rz(-2.7483181) q[0];
sx q[0];
rz(-1.9140592) q[0];
rz(-0.50499376) q[1];
sx q[1];
rz(-2.125449) q[1];
sx q[1];
rz(1.9379157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47984637) q[0];
sx q[0];
rz(-1.3239613) q[0];
sx q[0];
rz(1.4479475) q[0];
rz(0.92166047) q[2];
sx q[2];
rz(-1.1716915) q[2];
sx q[2];
rz(1.6726577) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.595049) q[1];
sx q[1];
rz(-2.1814934) q[1];
sx q[1];
rz(-2.2799019) q[1];
rz(-2.876128) q[3];
sx q[3];
rz(-2.0164312) q[3];
sx q[3];
rz(-1.3973941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0248854) q[2];
sx q[2];
rz(-1.025277) q[2];
sx q[2];
rz(-2.1825979) q[2];
rz(-0.2995019) q[3];
sx q[3];
rz(-0.12226573) q[3];
sx q[3];
rz(-1.449077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9129979) q[0];
sx q[0];
rz(-1.8916425) q[0];
sx q[0];
rz(2.6477497) q[0];
rz(0.85810581) q[1];
sx q[1];
rz(-1.162642) q[1];
sx q[1];
rz(1.3546622) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1345024) q[0];
sx q[0];
rz(-2.1190152) q[0];
sx q[0];
rz(-0.28315084) q[0];
x q[1];
rz(2.0628711) q[2];
sx q[2];
rz(-1.6920231) q[2];
sx q[2];
rz(-1.6693142) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.014381) q[1];
sx q[1];
rz(-0.76103044) q[1];
sx q[1];
rz(-3.0396023) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1936761) q[3];
sx q[3];
rz(-2.3661032) q[3];
sx q[3];
rz(0.38748565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4570423) q[2];
sx q[2];
rz(-2.5164618) q[2];
sx q[2];
rz(0.54614145) q[2];
rz(0.1977194) q[3];
sx q[3];
rz(-1.0309912) q[3];
sx q[3];
rz(2.8945967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2208743) q[0];
sx q[0];
rz(-1.9863167) q[0];
sx q[0];
rz(-2.0350631) q[0];
rz(-0.58440343) q[1];
sx q[1];
rz(-1.9677275) q[1];
sx q[1];
rz(2.1388785) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.913915) q[0];
sx q[0];
rz(-1.4105689) q[0];
sx q[0];
rz(-2.0100309) q[0];
x q[1];
rz(0.8876351) q[2];
sx q[2];
rz(-0.37345593) q[2];
sx q[2];
rz(-1.4591726) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.20368491) q[1];
sx q[1];
rz(-1.6069429) q[1];
sx q[1];
rz(-0.51236492) q[1];
rz(-1.9225305) q[3];
sx q[3];
rz(-2.1603843) q[3];
sx q[3];
rz(2.0812579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.57548412) q[2];
sx q[2];
rz(-2.2169952) q[2];
sx q[2];
rz(-1.3893348) q[2];
rz(0.53650457) q[3];
sx q[3];
rz(-1.6335231) q[3];
sx q[3];
rz(2.2947521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.093512) q[0];
sx q[0];
rz(-0.41963136) q[0];
sx q[0];
rz(-2.804948) q[0];
rz(3.0632784) q[1];
sx q[1];
rz(-1.9322194) q[1];
sx q[1];
rz(1.1557109) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5040249) q[0];
sx q[0];
rz(-2.2462566) q[0];
sx q[0];
rz(0.43677377) q[0];
rz(-pi) q[1];
rz(3.0750257) q[2];
sx q[2];
rz(-1.9835501) q[2];
sx q[2];
rz(-2.33605) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6558469) q[1];
sx q[1];
rz(-2.0585723) q[1];
sx q[1];
rz(2.9887245) q[1];
x q[2];
rz(-1.6011681) q[3];
sx q[3];
rz(-1.8874468) q[3];
sx q[3];
rz(2.3476072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.081850514) q[2];
sx q[2];
rz(-1.0826702) q[2];
sx q[2];
rz(-2.3183863) q[2];
rz(2.1246223) q[3];
sx q[3];
rz(-2.7401994) q[3];
sx q[3];
rz(-3.1034191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3391649) q[0];
sx q[0];
rz(-2.9186354) q[0];
sx q[0];
rz(-1.4270576) q[0];
rz(1.4786221) q[1];
sx q[1];
rz(-2.069811) q[1];
sx q[1];
rz(-0.85085416) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9768391) q[0];
sx q[0];
rz(-2.0434336) q[0];
sx q[0];
rz(-3.0125299) q[0];
x q[1];
rz(3.0548373) q[2];
sx q[2];
rz(-1.6977001) q[2];
sx q[2];
rz(-2.639304) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0180698) q[1];
sx q[1];
rz(-1.6241956) q[1];
sx q[1];
rz(0.52402871) q[1];
rz(-0.14031841) q[3];
sx q[3];
rz(-1.5068973) q[3];
sx q[3];
rz(1.7238703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6829546) q[2];
sx q[2];
rz(-1.7913603) q[2];
sx q[2];
rz(-0.82009912) q[2];
rz(-1.5470777) q[3];
sx q[3];
rz(-1.7087414) q[3];
sx q[3];
rz(1.1944176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9783258) q[0];
sx q[0];
rz(-1.5112725) q[0];
sx q[0];
rz(1.5164966) q[0];
rz(-0.77240472) q[1];
sx q[1];
rz(-2.518387) q[1];
sx q[1];
rz(1.0060681) q[1];
rz(-0.69140534) q[2];
sx q[2];
rz(-1.7476191) q[2];
sx q[2];
rz(-2.501001) q[2];
rz(-2.6670868) q[3];
sx q[3];
rz(-1.1474599) q[3];
sx q[3];
rz(-1.3808822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
