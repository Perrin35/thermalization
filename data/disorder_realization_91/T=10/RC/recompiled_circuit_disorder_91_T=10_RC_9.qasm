OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(5.7167238) q[0];
sx q[0];
rz(9.5958435) q[0];
rz(1.45362) q[1];
sx q[1];
rz(-0.34314081) q[1];
sx q[1];
rz(-1.3309825) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0199466) q[0];
sx q[0];
rz(-0.92187798) q[0];
sx q[0];
rz(-1.5050423) q[0];
rz(0.038161909) q[2];
sx q[2];
rz(-0.51617235) q[2];
sx q[2];
rz(-1.49522) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1316489) q[1];
sx q[1];
rz(-2.8899) q[1];
sx q[1];
rz(-1.5021055) q[1];
rz(-pi) q[2];
rz(-1.347723) q[3];
sx q[3];
rz(-2.2017456) q[3];
sx q[3];
rz(2.0032721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3657637) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(1.3295056) q[2];
rz(1.4154411) q[3];
sx q[3];
rz(-1.5751782) q[3];
sx q[3];
rz(-2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(2.1502894) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(2.9406722) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(2.8413049) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4703003) q[0];
sx q[0];
rz(-1.6912582) q[0];
sx q[0];
rz(-2.5380773) q[0];
rz(2.5710201) q[2];
sx q[2];
rz(-1.3010446) q[2];
sx q[2];
rz(0.74769339) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.092204658) q[1];
sx q[1];
rz(-0.37282473) q[1];
sx q[1];
rz(2.2224109) q[1];
rz(-pi) q[2];
rz(-1.3817915) q[3];
sx q[3];
rz(-1.4100037) q[3];
sx q[3];
rz(2.8652428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.369027) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(-1.0380113) q[2];
rz(2.299262) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(0.37030181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5623986) q[0];
sx q[0];
rz(-0.47214046) q[0];
sx q[0];
rz(-0.38811362) q[0];
rz(0.072470486) q[1];
sx q[1];
rz(-1.7150755) q[1];
sx q[1];
rz(2.8252576) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2280699) q[0];
sx q[0];
rz(-1.5689578) q[0];
sx q[0];
rz(-1.3257922) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6109654) q[2];
sx q[2];
rz(-0.66662153) q[2];
sx q[2];
rz(-1.9463469) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2410779) q[1];
sx q[1];
rz(-1.1763651) q[1];
sx q[1];
rz(-2.6532252) q[1];
rz(-pi) q[2];
rz(-2.2736336) q[3];
sx q[3];
rz(-0.43324019) q[3];
sx q[3];
rz(0.058335282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1091653) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(-0.16080984) q[2];
rz(-3.0155449) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(1.0236615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2407103) q[0];
sx q[0];
rz(-0.70403376) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(1.7968934) q[1];
sx q[1];
rz(-0.99935499) q[1];
sx q[1];
rz(1.5136738) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31878372) q[0];
sx q[0];
rz(-1.6167287) q[0];
sx q[0];
rz(-2.4687693) q[0];
rz(1.4430181) q[2];
sx q[2];
rz(-0.56782349) q[2];
sx q[2];
rz(2.36433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0878144) q[1];
sx q[1];
rz(-2.6930241) q[1];
sx q[1];
rz(1.5320369) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73174814) q[3];
sx q[3];
rz(-1.8879315) q[3];
sx q[3];
rz(-1.5665311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7164798) q[2];
sx q[2];
rz(-2.1322865) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(2.7246357) q[3];
sx q[3];
rz(-2.0988393) q[3];
sx q[3];
rz(0.52880374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.058218) q[0];
sx q[0];
rz(-0.26979065) q[0];
sx q[0];
rz(1.8378687) q[0];
rz(2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(0.051503332) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6131825) q[0];
sx q[0];
rz(-0.15325704) q[0];
sx q[0];
rz(-0.86092237) q[0];
x q[1];
rz(1.0837469) q[2];
sx q[2];
rz(-0.32099989) q[2];
sx q[2];
rz(0.1333065) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0609378) q[1];
sx q[1];
rz(-0.47532156) q[1];
sx q[1];
rz(1.2259543) q[1];
rz(-pi) q[2];
rz(-3.0052161) q[3];
sx q[3];
rz(-0.73835056) q[3];
sx q[3];
rz(-0.34463681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(0.55001843) q[3];
sx q[3];
rz(-0.23967448) q[3];
sx q[3];
rz(1.1365183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0710058) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(0.20198527) q[0];
rz(0.96549353) q[1];
sx q[1];
rz(-2.1172724) q[1];
sx q[1];
rz(3.0583256) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4250454) q[0];
sx q[0];
rz(-1.7444567) q[0];
sx q[0];
rz(-2.9437149) q[0];
x q[1];
rz(-1.8385356) q[2];
sx q[2];
rz(-1.8362852) q[2];
sx q[2];
rz(0.43648411) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80379936) q[1];
sx q[1];
rz(-2.7593136) q[1];
sx q[1];
rz(1.2804968) q[1];
x q[2];
rz(1.8694599) q[3];
sx q[3];
rz(-1.0335869) q[3];
sx q[3];
rz(2.1512973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.10963708) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(1.7588245) q[2];
rz(0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4270585) q[0];
sx q[0];
rz(-0.80474168) q[0];
sx q[0];
rz(-0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(-2.0239963) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2876802) q[0];
sx q[0];
rz(-0.84830647) q[0];
sx q[0];
rz(1.393909) q[0];
x q[1];
rz(2.2824077) q[2];
sx q[2];
rz(-1.7297598) q[2];
sx q[2];
rz(0.32260103) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0301789) q[1];
sx q[1];
rz(-1.7708659) q[1];
sx q[1];
rz(-2.5740037) q[1];
rz(-2.6071129) q[3];
sx q[3];
rz(-1.3325053) q[3];
sx q[3];
rz(2.0508545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.08427944) q[2];
sx q[2];
rz(-1.7417615) q[2];
sx q[2];
rz(-2.9220667) q[2];
rz(-2.7548742) q[3];
sx q[3];
rz(-2.3733449) q[3];
sx q[3];
rz(0.99532551) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052208386) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(-0.21729939) q[0];
rz(3.1106588) q[1];
sx q[1];
rz(-0.63806454) q[1];
sx q[1];
rz(2.4826179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156408) q[0];
sx q[0];
rz(-0.77227393) q[0];
sx q[0];
rz(3.1053567) q[0];
x q[1];
rz(0.73924139) q[2];
sx q[2];
rz(-2.2695702) q[2];
sx q[2];
rz(3.0657363) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1606071) q[1];
sx q[1];
rz(-1.7874582) q[1];
sx q[1];
rz(0.76088455) q[1];
rz(-pi) q[2];
rz(-0.92774763) q[3];
sx q[3];
rz(-1.5156931) q[3];
sx q[3];
rz(0.4004713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7028246) q[2];
sx q[2];
rz(-1.7610901) q[2];
sx q[2];
rz(-2.8213815) q[2];
rz(1.5252339) q[3];
sx q[3];
rz(-1.1956918) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3808688) q[0];
sx q[0];
rz(-0.18146935) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(-1.0702417) q[1];
sx q[1];
rz(-1.8990592) q[1];
sx q[1];
rz(1.9314996) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4142128) q[0];
sx q[0];
rz(-0.41398898) q[0];
sx q[0];
rz(-1.0340286) q[0];
rz(-0.25173431) q[2];
sx q[2];
rz(-2.1605957) q[2];
sx q[2];
rz(-0.41433197) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6283419) q[1];
sx q[1];
rz(-1.8645789) q[1];
sx q[1];
rz(-0.88900868) q[1];
rz(2.6415778) q[3];
sx q[3];
rz(-2.4588636) q[3];
sx q[3];
rz(2.6422215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.575763) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(0.56662095) q[2];
rz(-0.9295272) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41314769) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(1.43191) q[0];
rz(-2.6152949) q[1];
sx q[1];
rz(-2.6142575) q[1];
sx q[1];
rz(0.79968232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0678588) q[0];
sx q[0];
rz(-1.6556544) q[0];
sx q[0];
rz(-1.2220864) q[0];
rz(-pi) q[1];
rz(1.6522371) q[2];
sx q[2];
rz(-1.4105721) q[2];
sx q[2];
rz(2.3022431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4754776) q[1];
sx q[1];
rz(-2.0787422) q[1];
sx q[1];
rz(2.75027) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6639789) q[3];
sx q[3];
rz(-2.9778746) q[3];
sx q[3];
rz(1.175566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2852823) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(-0.69520673) q[2];
rz(0.55784145) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0859062) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(-3.1391115) q[2];
sx q[2];
rz(-0.33591349) q[2];
sx q[2];
rz(-2.9813067) q[2];
rz(-1.7351032) q[3];
sx q[3];
rz(-0.77948979) q[3];
sx q[3];
rz(-0.43850552) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
