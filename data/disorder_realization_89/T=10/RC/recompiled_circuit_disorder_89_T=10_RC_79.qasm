OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(0.29456079) q[0];
rz(-2.8566868) q[1];
sx q[1];
rz(-2.6309738) q[1];
sx q[1];
rz(2.7198305) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.598782) q[0];
sx q[0];
rz(-3.071975) q[0];
sx q[0];
rz(1.1845864) q[0];
rz(-pi) q[1];
rz(-1.6099036) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(-1.4456911) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0114853) q[1];
sx q[1];
rz(-1.3890146) q[1];
sx q[1];
rz(2.6891522) q[1];
rz(-pi) q[2];
rz(-2.9362039) q[3];
sx q[3];
rz(-2.5377512) q[3];
sx q[3];
rz(0.83128923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1224147) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(-3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(2.2825867) q[0];
rz(2.7711218) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252936) q[0];
sx q[0];
rz(-1.787961) q[0];
sx q[0];
rz(-0.3152245) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53441647) q[2];
sx q[2];
rz(-2.538531) q[2];
sx q[2];
rz(-0.39819983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6795923) q[1];
sx q[1];
rz(-0.781171) q[1];
sx q[1];
rz(2.6189234) q[1];
rz(1.8584077) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(0.19954296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(2.583288) q[2];
rz(-2.1022508) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(2.6168183) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(-2.4257461) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(-1.3735501) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4931902) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(-0.52669749) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85767304) q[2];
sx q[2];
rz(-1.8048986) q[2];
sx q[2];
rz(0.26299325) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9436685) q[1];
sx q[1];
rz(-2.1197882) q[1];
sx q[1];
rz(-0.82358574) q[1];
rz(-1.8009637) q[3];
sx q[3];
rz(-1.3535415) q[3];
sx q[3];
rz(1.0297071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0107161) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(1.019657) q[2];
rz(-2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(-2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402128) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(-2.3024978) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2903039) q[0];
sx q[0];
rz(-1.5595058) q[0];
sx q[0];
rz(-1.4615061) q[0];
x q[1];
rz(1.9539815) q[2];
sx q[2];
rz(-1.9028579) q[2];
sx q[2];
rz(3.1021169) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0415045) q[1];
sx q[1];
rz(-1.8309635) q[1];
sx q[1];
rz(2.2799904) q[1];
rz(-2.2545492) q[3];
sx q[3];
rz(-2.4409557) q[3];
sx q[3];
rz(-2.5168583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6976167) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(2.9042517) q[2];
rz(0.90732968) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(-1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6326555) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.5326112) q[0];
rz(-2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(-1.8283432) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35996138) q[0];
sx q[0];
rz(-0.93802035) q[0];
sx q[0];
rz(-2.7347793) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1336742) q[2];
sx q[2];
rz(-0.067069947) q[2];
sx q[2];
rz(-1.1941393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.207706) q[1];
sx q[1];
rz(-1.943214) q[1];
sx q[1];
rz(-1.1697342) q[1];
rz(-pi) q[2];
rz(-0.58532183) q[3];
sx q[3];
rz(-0.96671852) q[3];
sx q[3];
rz(2.4620591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0017172) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(-0.66429794) q[2];
rz(0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(-3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649081) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(3.0867807) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(-0.48318133) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9478215) q[0];
sx q[0];
rz(-1.3894488) q[0];
sx q[0];
rz(3.0788455) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1088164) q[2];
sx q[2];
rz(-2.0332697) q[2];
sx q[2];
rz(1.0810766) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.74181611) q[1];
sx q[1];
rz(-0.92405926) q[1];
sx q[1];
rz(-2.497655) q[1];
rz(-pi) q[2];
rz(3.0275434) q[3];
sx q[3];
rz(-1.9012791) q[3];
sx q[3];
rz(1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.15423916) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(-3.0563291) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(0.15032642) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(-2.1869587) q[0];
rz(-0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(0.25209299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1958256) q[0];
sx q[0];
rz(-2.6951163) q[0];
sx q[0];
rz(2.8004942) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3605395) q[2];
sx q[2];
rz(-2.1456246) q[2];
sx q[2];
rz(0.99036723) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.39602173) q[1];
sx q[1];
rz(-0.77116291) q[1];
sx q[1];
rz(2.1019756) q[1];
rz(2.2456456) q[3];
sx q[3];
rz(-0.5657256) q[3];
sx q[3];
rz(-1.56324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(0.90448109) q[2];
rz(-2.1379743) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(-1.4138387) q[0];
rz(-0.66043234) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(1.3716912) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8577514) q[0];
sx q[0];
rz(-0.62566602) q[0];
sx q[0];
rz(1.4197423) q[0];
x q[1];
rz(-2.9210864) q[2];
sx q[2];
rz(-2.3683511) q[2];
sx q[2];
rz(-0.79794937) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.614537) q[1];
sx q[1];
rz(-1.0053047) q[1];
sx q[1];
rz(-1.2177699) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7226108) q[3];
sx q[3];
rz(-1.2137128) q[3];
sx q[3];
rz(2.9861772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-0.39789847) q[2];
rz(-0.49063101) q[3];
sx q[3];
rz(-1.5964973) q[3];
sx q[3];
rz(1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21815498) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(0.27059069) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9070248) q[0];
sx q[0];
rz(-1.3493291) q[0];
sx q[0];
rz(1.7937167) q[0];
rz(-1.2241237) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(-0.62435645) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8671659) q[1];
sx q[1];
rz(-1.567489) q[1];
sx q[1];
rz(-0.13722384) q[1];
rz(-1.0911646) q[3];
sx q[3];
rz(-1.661146) q[3];
sx q[3];
rz(-2.3144212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33971912) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(1.489893) q[2];
rz(-2.4387032) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-0.15429601) q[0];
rz(0.94296304) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(2.399209) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48001227) q[0];
sx q[0];
rz(-0.98012692) q[0];
sx q[0];
rz(2.3175879) q[0];
rz(-pi) q[1];
rz(-0.97147271) q[2];
sx q[2];
rz(-0.68442548) q[2];
sx q[2];
rz(1.7873834) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2652475) q[1];
sx q[1];
rz(-0.34595385) q[1];
sx q[1];
rz(-0.032851263) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14245716) q[3];
sx q[3];
rz(-1.5995306) q[3];
sx q[3];
rz(1.8491668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2877038) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29006526) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(-1.3416946) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(-2.3201597) q[2];
sx q[2];
rz(-1.8100304) q[2];
sx q[2];
rz(2.9678154) q[2];
rz(-0.7704173) q[3];
sx q[3];
rz(-2.9478248) q[3];
sx q[3];
rz(-0.41268681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
