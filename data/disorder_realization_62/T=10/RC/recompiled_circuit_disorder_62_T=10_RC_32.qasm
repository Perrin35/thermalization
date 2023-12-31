OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(-2.9641889) q[0];
sx q[0];
rz(2.0071964) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(-2.1048574) q[1];
sx q[1];
rz(-0.66361767) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3463319) q[0];
sx q[0];
rz(-1.5471317) q[0];
sx q[0];
rz(-2.6803826) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3637949) q[2];
sx q[2];
rz(-1.8282837) q[2];
sx q[2];
rz(2.4553026) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.2634969) q[1];
sx q[1];
rz(-2.7090008) q[1];
sx q[1];
rz(0.90579512) q[1];
x q[2];
rz(1.8144242) q[3];
sx q[3];
rz(-1.4654136) q[3];
sx q[3];
rz(1.5308612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.87876451) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(-0.051068548) q[2];
rz(2.5845394) q[3];
sx q[3];
rz(-0.80018187) q[3];
sx q[3];
rz(-1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(0.59659514) q[0];
rz(-2.3157628) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.2260431) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37218371) q[0];
sx q[0];
rz(-2.6499977) q[0];
sx q[0];
rz(-0.43222897) q[0];
rz(-pi) q[1];
rz(-2.3184899) q[2];
sx q[2];
rz(-0.96379333) q[2];
sx q[2];
rz(-2.5603103) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5920168) q[1];
sx q[1];
rz(-1.5703652) q[1];
sx q[1];
rz(-1.3509343) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1434104) q[3];
sx q[3];
rz(-1.5628353) q[3];
sx q[3];
rz(1.5496467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3423959) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(-2.8105695) q[2];
rz(0.80667574) q[3];
sx q[3];
rz(-1.9162063) q[3];
sx q[3];
rz(1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.230343) q[0];
sx q[0];
rz(1.249041) q[0];
rz(-3.0535835) q[1];
sx q[1];
rz(-1.1227337) q[1];
sx q[1];
rz(-1.0294611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3345966) q[0];
sx q[0];
rz(-2.1801729) q[0];
sx q[0];
rz(0.2683123) q[0];
rz(2.4918409) q[2];
sx q[2];
rz(-1.3403112) q[2];
sx q[2];
rz(0.31313716) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6368235) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(1.7117281) q[1];
rz(-pi) q[2];
rz(2.401628) q[3];
sx q[3];
rz(-1.660534) q[3];
sx q[3];
rz(-1.2147853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.133698) q[2];
sx q[2];
rz(-1.4174613) q[2];
sx q[2];
rz(0.50764817) q[2];
rz(-1.3890022) q[3];
sx q[3];
rz(-2.8116083) q[3];
sx q[3];
rz(-1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65748173) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(1.5456276) q[0];
rz(2.0987299) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(-1.5159336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.078331) q[0];
sx q[0];
rz(-2.445979) q[0];
sx q[0];
rz(-0.22380933) q[0];
x q[1];
rz(-1.707294) q[2];
sx q[2];
rz(-2.578306) q[2];
sx q[2];
rz(-2.0493281) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8856636) q[1];
sx q[1];
rz(-2.4293578) q[1];
sx q[1];
rz(0.42339143) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5528615) q[3];
sx q[3];
rz(-1.3874467) q[3];
sx q[3];
rz(-0.4962894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3778014) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-2.5989792) q[3];
sx q[3];
rz(2.1550089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7871053) q[0];
sx q[0];
rz(-1.1802477) q[0];
sx q[0];
rz(-1.6217344) q[0];
rz(2.5095818) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(0.89486665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2762404) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(-2.4549237) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.61770265) q[2];
sx q[2];
rz(-2.3139944) q[2];
sx q[2];
rz(-0.4862116) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8708961) q[1];
sx q[1];
rz(-1.1333229) q[1];
sx q[1];
rz(0.84769627) q[1];
x q[2];
rz(-1.5630866) q[3];
sx q[3];
rz(-1.8572516) q[3];
sx q[3];
rz(2.8849998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.616509) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(1.1425225) q[2];
rz(-2.3948495) q[3];
sx q[3];
rz(-1.8871504) q[3];
sx q[3];
rz(-1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557945) q[0];
sx q[0];
rz(-0.32145158) q[0];
sx q[0];
rz(-1.7640132) q[0];
rz(0.47239834) q[1];
sx q[1];
rz(-0.51858416) q[1];
sx q[1];
rz(-0.46498743) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0156292) q[0];
sx q[0];
rz(-2.0032126) q[0];
sx q[0];
rz(-2.7727491) q[0];
rz(-pi) q[1];
rz(1.497252) q[2];
sx q[2];
rz(-1.911517) q[2];
sx q[2];
rz(-1.3982915) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7558414) q[1];
sx q[1];
rz(-0.63703905) q[1];
sx q[1];
rz(-1.6511276) q[1];
rz(-pi) q[2];
rz(2.898071) q[3];
sx q[3];
rz(-0.72391073) q[3];
sx q[3];
rz(0.078725423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.650699) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(2.6965551) q[2];
rz(2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(0.26708189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0141107) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.7215464) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-0.61444608) q[1];
sx q[1];
rz(2.9856317) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6316846) q[0];
sx q[0];
rz(-1.317306) q[0];
sx q[0];
rz(1.3939199) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1440802) q[2];
sx q[2];
rz(-2.0442171) q[2];
sx q[2];
rz(2.0063426) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5739463) q[1];
sx q[1];
rz(-1.5168377) q[1];
sx q[1];
rz(1.3957363) q[1];
rz(-pi) q[2];
rz(2.3994) q[3];
sx q[3];
rz(-2.4626197) q[3];
sx q[3];
rz(1.9051444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.069313958) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(2.8159451) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(1.5163039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9496562) q[0];
sx q[0];
rz(-0.40238109) q[0];
sx q[0];
rz(-0.33777133) q[0];
rz(1.0900963) q[1];
sx q[1];
rz(-2.1665159) q[1];
sx q[1];
rz(-0.24857323) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2598495) q[0];
sx q[0];
rz(-1.8200841) q[0];
sx q[0];
rz(-2.6134406) q[0];
x q[1];
rz(0.48034251) q[2];
sx q[2];
rz(-2.0311653) q[2];
sx q[2];
rz(0.72887052) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.92330248) q[1];
sx q[1];
rz(-0.81070886) q[1];
sx q[1];
rz(-0.87358012) q[1];
x q[2];
rz(1.400984) q[3];
sx q[3];
rz(-1.4949993) q[3];
sx q[3];
rz(0.71789391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2934072) q[2];
sx q[2];
rz(-2.6082787) q[2];
sx q[2];
rz(0.08671134) q[2];
rz(0.48197204) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6329353) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(-0.28276643) q[0];
rz(2.4400318) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(-1.823002) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628172) q[0];
sx q[0];
rz(-1.6512198) q[0];
sx q[0];
rz(1.1286939) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65638541) q[2];
sx q[2];
rz(-1.0644541) q[2];
sx q[2];
rz(2.2040747) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5732167) q[1];
sx q[1];
rz(-1.425256) q[1];
sx q[1];
rz(2.2578866) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1845469) q[3];
sx q[3];
rz(-0.85856122) q[3];
sx q[3];
rz(-1.4854747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(0.39917699) q[2];
rz(2.2579851) q[3];
sx q[3];
rz(-1.4727605) q[3];
sx q[3];
rz(1.9201027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3783962) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(-1.5989074) q[0];
rz(1.0653161) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(-1.8803966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1051837) q[0];
sx q[0];
rz(-1.1622218) q[0];
sx q[0];
rz(1.6856134) q[0];
x q[1];
rz(-0.59826675) q[2];
sx q[2];
rz(-1.5891799) q[2];
sx q[2];
rz(-0.53668864) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.48325) q[1];
sx q[1];
rz(-1.5066506) q[1];
sx q[1];
rz(2.642753) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4576549) q[3];
sx q[3];
rz(-1.4468907) q[3];
sx q[3];
rz(-0.77139664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5836872) q[2];
sx q[2];
rz(-0.62290278) q[2];
sx q[2];
rz(2.5718001) q[2];
rz(1.9231046) q[3];
sx q[3];
rz(-0.64703882) q[3];
sx q[3];
rz(-0.56263721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31310836) q[0];
sx q[0];
rz(-0.88043558) q[0];
sx q[0];
rz(1.2783929) q[0];
rz(0.60824153) q[1];
sx q[1];
rz(-2.6651762) q[1];
sx q[1];
rz(-2.6574635) q[1];
rz(-2.1424978) q[2];
sx q[2];
rz(-1.5928762) q[2];
sx q[2];
rz(1.5230509) q[2];
rz(-0.20523397) q[3];
sx q[3];
rz(-0.60094613) q[3];
sx q[3];
rz(-0.080106674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
