OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3818504) q[0];
sx q[0];
rz(-0.83431017) q[0];
sx q[0];
rz(-0.068339737) q[0];
rz(-0.66032687) q[1];
sx q[1];
rz(-0.84815174) q[1];
sx q[1];
rz(-3.1037722) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4662287) q[0];
sx q[0];
rz(-2.9328406) q[0];
sx q[0];
rz(-1.1842313) q[0];
x q[1];
rz(-2.2322466) q[2];
sx q[2];
rz(-1.6644018) q[2];
sx q[2];
rz(2.0441165) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6112069) q[1];
sx q[1];
rz(-1.5198277) q[1];
sx q[1];
rz(0.74688046) q[1];
rz(-1.9170403) q[3];
sx q[3];
rz(-1.4284705) q[3];
sx q[3];
rz(-3.0179253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7906856) q[2];
sx q[2];
rz(-0.51351341) q[2];
sx q[2];
rz(1.2834056) q[2];
rz(-3.0154199) q[3];
sx q[3];
rz(-1.4161371) q[3];
sx q[3];
rz(-3.0509907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44777563) q[0];
sx q[0];
rz(-0.87647804) q[0];
sx q[0];
rz(-0.59666657) q[0];
rz(-1.5555351) q[1];
sx q[1];
rz(-1.3417473) q[1];
sx q[1];
rz(-1.3635051) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039149337) q[0];
sx q[0];
rz(-1.9958695) q[0];
sx q[0];
rz(-2.5545679) q[0];
x q[1];
rz(-3.0599669) q[2];
sx q[2];
rz(-0.55742369) q[2];
sx q[2];
rz(2.4068085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5433568) q[1];
sx q[1];
rz(-2.1992866) q[1];
sx q[1];
rz(-1.9176107) q[1];
rz(0.73804654) q[3];
sx q[3];
rz(-0.9934721) q[3];
sx q[3];
rz(-0.20418024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3229225) q[2];
sx q[2];
rz(-2.1408036) q[2];
sx q[2];
rz(-0.73454109) q[2];
rz(-2.5143886) q[3];
sx q[3];
rz(-1.5137129) q[3];
sx q[3];
rz(2.396615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4194141) q[0];
sx q[0];
rz(-2.3086771) q[0];
sx q[0];
rz(-2.869379) q[0];
rz(0.84725562) q[1];
sx q[1];
rz(-1.3191185) q[1];
sx q[1];
rz(-0.31262696) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6486559) q[0];
sx q[0];
rz(-1.3515633) q[0];
sx q[0];
rz(-1.6256902) q[0];
rz(-pi) q[1];
rz(2.0414511) q[2];
sx q[2];
rz(-1.0725642) q[2];
sx q[2];
rz(2.2217896) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9885013) q[1];
sx q[1];
rz(-1.3107436) q[1];
sx q[1];
rz(-0.071916332) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7496495) q[3];
sx q[3];
rz(-2.3392945) q[3];
sx q[3];
rz(0.18303579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0212038) q[2];
sx q[2];
rz(-2.618232) q[2];
sx q[2];
rz(0.81494251) q[2];
rz(-0.96873823) q[3];
sx q[3];
rz(-1.0453984) q[3];
sx q[3];
rz(2.708784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.3880436) q[0];
sx q[0];
rz(-0.20064813) q[0];
sx q[0];
rz(2.5182305) q[0];
rz(0.81758824) q[1];
sx q[1];
rz(-0.31660429) q[1];
sx q[1];
rz(-0.049291074) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3885956) q[0];
sx q[0];
rz(-1.5224644) q[0];
sx q[0];
rz(-2.3301793) q[0];
rz(0.85929112) q[2];
sx q[2];
rz(-1.8154732) q[2];
sx q[2];
rz(-1.7161075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4781487) q[1];
sx q[1];
rz(-1.2423007) q[1];
sx q[1];
rz(-2.8847787) q[1];
rz(2.2634099) q[3];
sx q[3];
rz(-2.4947824) q[3];
sx q[3];
rz(1.1991024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7307044) q[2];
sx q[2];
rz(-1.5894019) q[2];
sx q[2];
rz(-0.98497406) q[2];
rz(-2.2385521) q[3];
sx q[3];
rz(-1.1629546) q[3];
sx q[3];
rz(0.27339098) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34489283) q[0];
sx q[0];
rz(-1.5364237) q[0];
sx q[0];
rz(3.0539736) q[0];
rz(2.9852729) q[1];
sx q[1];
rz(-2.5904398) q[1];
sx q[1];
rz(2.2706251) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6773416) q[0];
sx q[0];
rz(-0.28872492) q[0];
sx q[0];
rz(-3.0407453) q[0];
rz(1.1768612) q[2];
sx q[2];
rz(-1.7125868) q[2];
sx q[2];
rz(2.1243492) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85775447) q[1];
sx q[1];
rz(-0.72099599) q[1];
sx q[1];
rz(-0.43197076) q[1];
rz(-pi) q[2];
rz(2.3351339) q[3];
sx q[3];
rz(-1.4440923) q[3];
sx q[3];
rz(-2.5106631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.093420371) q[2];
sx q[2];
rz(-2.293736) q[2];
sx q[2];
rz(-0.35287738) q[2];
rz(2.2757163) q[3];
sx q[3];
rz(-0.70169774) q[3];
sx q[3];
rz(-0.29233366) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6784994) q[0];
sx q[0];
rz(-1.1266288) q[0];
sx q[0];
rz(-3.034806) q[0];
rz(1.9550025) q[1];
sx q[1];
rz(-1.0266961) q[1];
sx q[1];
rz(-1.0245163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9232641) q[0];
sx q[0];
rz(-0.54696333) q[0];
sx q[0];
rz(2.5640019) q[0];
x q[1];
rz(2.5249135) q[2];
sx q[2];
rz(-2.4384535) q[2];
sx q[2];
rz(2.6577735) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.34895996) q[1];
sx q[1];
rz(-1.0981202) q[1];
sx q[1];
rz(2.2036168) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.939219) q[3];
sx q[3];
rz(-2.0344067) q[3];
sx q[3];
rz(1.1045477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7496877) q[2];
sx q[2];
rz(-1.2140032) q[2];
sx q[2];
rz(2.270703) q[2];
rz(-1.8093367) q[3];
sx q[3];
rz(-2.199316) q[3];
sx q[3];
rz(1.4107305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9695327) q[0];
sx q[0];
rz(-1.7344069) q[0];
sx q[0];
rz(2.5906738) q[0];
rz(-2.6761966) q[1];
sx q[1];
rz(-2.4702563) q[1];
sx q[1];
rz(0.23682061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.119334) q[0];
sx q[0];
rz(-0.10611457) q[0];
sx q[0];
rz(-1.2073713) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15010712) q[2];
sx q[2];
rz(-1.7068212) q[2];
sx q[2];
rz(1.0667691) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.782978) q[1];
sx q[1];
rz(-1.8404984) q[1];
sx q[1];
rz(3.0637834) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3051466) q[3];
sx q[3];
rz(-1.340938) q[3];
sx q[3];
rz(3.0695855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44234309) q[2];
sx q[2];
rz(-0.68168679) q[2];
sx q[2];
rz(-2.701475) q[2];
rz(-1.0951428) q[3];
sx q[3];
rz(-1.6710072) q[3];
sx q[3];
rz(1.7949036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(2.8286164) q[0];
sx q[0];
rz(-1.3422817) q[0];
sx q[0];
rz(3.112088) q[0];
rz(-0.94738952) q[1];
sx q[1];
rz(-2.3209929) q[1];
sx q[1];
rz(3.0227919) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4122075) q[0];
sx q[0];
rz(-0.31353024) q[0];
sx q[0];
rz(0.14051147) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88840862) q[2];
sx q[2];
rz(-2.4783122) q[2];
sx q[2];
rz(2.1036489) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48092948) q[1];
sx q[1];
rz(-2.5813563) q[1];
sx q[1];
rz(-2.5653097) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8700065) q[3];
sx q[3];
rz(-0.51338235) q[3];
sx q[3];
rz(-1.3083003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7416731) q[2];
sx q[2];
rz(-0.46367773) q[2];
sx q[2];
rz(-1.5681533) q[2];
rz(-1.167477) q[3];
sx q[3];
rz(-1.7220595) q[3];
sx q[3];
rz(1.2742111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2385999) q[0];
sx q[0];
rz(-0.82505834) q[0];
sx q[0];
rz(0.15326823) q[0];
rz(2.080147) q[1];
sx q[1];
rz(-0.57466424) q[1];
sx q[1];
rz(-1.1538039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.783219) q[0];
sx q[0];
rz(-1.3010539) q[0];
sx q[0];
rz(-0.37958919) q[0];
rz(-pi) q[1];
rz(0.93038721) q[2];
sx q[2];
rz(-0.79158917) q[2];
sx q[2];
rz(1.1861578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.02281636) q[1];
sx q[1];
rz(-2.9252242) q[1];
sx q[1];
rz(-1.3092036) q[1];
rz(-pi) q[2];
rz(0.88345248) q[3];
sx q[3];
rz(-1.4769819) q[3];
sx q[3];
rz(2.1283619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8447421) q[2];
sx q[2];
rz(-1.1721609) q[2];
sx q[2];
rz(1.6142169) q[2];
rz(-0.99689233) q[3];
sx q[3];
rz(-1.4930054) q[3];
sx q[3];
rz(-0.87866384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.24755724) q[0];
sx q[0];
rz(-1.0972801) q[0];
sx q[0];
rz(2.9515008) q[0];
rz(2.4709573) q[1];
sx q[1];
rz(-1.1963444) q[1];
sx q[1];
rz(-1.488283) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22083902) q[0];
sx q[0];
rz(-1.5399884) q[0];
sx q[0];
rz(3.0160745) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52942099) q[2];
sx q[2];
rz(-0.501907) q[2];
sx q[2];
rz(-1.5159964) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6675122) q[1];
sx q[1];
rz(-1.0682567) q[1];
sx q[1];
rz(-2.4221592) q[1];
rz(-pi) q[2];
rz(-0.87741239) q[3];
sx q[3];
rz(-1.1023695) q[3];
sx q[3];
rz(-1.7978158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1378479) q[2];
sx q[2];
rz(-0.90128428) q[2];
sx q[2];
rz(-0.89912644) q[2];
rz(-2.6265465) q[3];
sx q[3];
rz(-0.47596541) q[3];
sx q[3];
rz(2.9639444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0638194) q[0];
sx q[0];
rz(-1.3615006) q[0];
sx q[0];
rz(2.5831945) q[0];
rz(0.28941119) q[1];
sx q[1];
rz(-0.90129539) q[1];
sx q[1];
rz(1.7064066) q[1];
rz(0.70984798) q[2];
sx q[2];
rz(-1.3550497) q[2];
sx q[2];
rz(1.4646127) q[2];
rz(-0.12712052) q[3];
sx q[3];
rz(-1.8416497) q[3];
sx q[3];
rz(1.9491378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
