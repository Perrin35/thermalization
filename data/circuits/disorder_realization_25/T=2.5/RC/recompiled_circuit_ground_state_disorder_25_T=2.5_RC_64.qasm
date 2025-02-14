OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.175617) q[0];
sx q[0];
rz(4.7369851) q[0];
sx q[0];
rz(10.935258) q[0];
rz(1.0881967) q[1];
sx q[1];
rz(-1.3047941) q[1];
sx q[1];
rz(2.3545177) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6841976) q[0];
sx q[0];
rz(-1.6826864) q[0];
sx q[0];
rz(0.68023139) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6252615) q[2];
sx q[2];
rz(-1.174841) q[2];
sx q[2];
rz(1.2690282) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0144498) q[1];
sx q[1];
rz(-2.9269955) q[1];
sx q[1];
rz(-0.17226528) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3863411) q[3];
sx q[3];
rz(-0.77338615) q[3];
sx q[3];
rz(0.36808792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2934908) q[2];
sx q[2];
rz(-3.0300671) q[2];
sx q[2];
rz(2.9254986) q[2];
rz(0.50348336) q[3];
sx q[3];
rz(-1.2516021) q[3];
sx q[3];
rz(0.066472806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1450495) q[0];
sx q[0];
rz(-1.1798877) q[0];
sx q[0];
rz(-2.766372) q[0];
rz(0.61620617) q[1];
sx q[1];
rz(-1.0169949) q[1];
sx q[1];
rz(0.18888758) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66697648) q[0];
sx q[0];
rz(-0.89578264) q[0];
sx q[0];
rz(2.6029909) q[0];
rz(-pi) q[1];
rz(-2.3023841) q[2];
sx q[2];
rz(-1.1611934) q[2];
sx q[2];
rz(1.6466779) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1950043) q[1];
sx q[1];
rz(-0.71041162) q[1];
sx q[1];
rz(-1.2571774) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6140574) q[3];
sx q[3];
rz(-1.9200824) q[3];
sx q[3];
rz(-2.7852102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6836493) q[2];
sx q[2];
rz(-1.8234437) q[2];
sx q[2];
rz(1.9981492) q[2];
rz(-2.5055366) q[3];
sx q[3];
rz(-0.053074107) q[3];
sx q[3];
rz(1.5422356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3627593) q[0];
sx q[0];
rz(-2.9017359) q[0];
sx q[0];
rz(-1.1676316) q[0];
rz(-3.0150343) q[1];
sx q[1];
rz(-1.5153706) q[1];
sx q[1];
rz(0.57356858) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99988692) q[0];
sx q[0];
rz(-0.83214271) q[0];
sx q[0];
rz(-0.12267273) q[0];
rz(0.36678183) q[2];
sx q[2];
rz(-1.8213118) q[2];
sx q[2];
rz(2.9652852) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6918317) q[1];
sx q[1];
rz(-2.0013325) q[1];
sx q[1];
rz(-0.052816347) q[1];
x q[2];
rz(2.4819991) q[3];
sx q[3];
rz(-1.9232188) q[3];
sx q[3];
rz(-3.0291758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48587376) q[2];
sx q[2];
rz(-1.8334917) q[2];
sx q[2];
rz(1.5029079) q[2];
rz(-1.8466628) q[3];
sx q[3];
rz(-0.27532268) q[3];
sx q[3];
rz(0.82726971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5824222) q[0];
sx q[0];
rz(-0.14635135) q[0];
sx q[0];
rz(0.91648066) q[0];
rz(0.16042635) q[1];
sx q[1];
rz(-1.285099) q[1];
sx q[1];
rz(-0.697335) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19577414) q[0];
sx q[0];
rz(-1.8911288) q[0];
sx q[0];
rz(2.8813502) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.38390145) q[2];
sx q[2];
rz(-2.4426201) q[2];
sx q[2];
rz(2.0052103) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26665822) q[1];
sx q[1];
rz(-1.8968079) q[1];
sx q[1];
rz(0.74649109) q[1];
x q[2];
rz(3.0352108) q[3];
sx q[3];
rz(-0.13088317) q[3];
sx q[3];
rz(1.037055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15630284) q[2];
sx q[2];
rz(-2.6142945) q[2];
sx q[2];
rz(2.018833) q[2];
rz(1.5491693) q[3];
sx q[3];
rz(-1.6808108) q[3];
sx q[3];
rz(2.7482225) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8321946) q[0];
sx q[0];
rz(-3.0351312) q[0];
sx q[0];
rz(1.1291809) q[0];
rz(0.87942266) q[1];
sx q[1];
rz(-1.3112473) q[1];
sx q[1];
rz(2.5130491) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1482233) q[0];
sx q[0];
rz(-1.4472479) q[0];
sx q[0];
rz(-0.11798239) q[0];
rz(-pi) q[1];
rz(1.0209924) q[2];
sx q[2];
rz(-0.98787427) q[2];
sx q[2];
rz(-1.1958808) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.4820833) q[1];
sx q[1];
rz(-1.0190367) q[1];
sx q[1];
rz(3.0943046) q[1];
rz(-1.699963) q[3];
sx q[3];
rz(-0.74366513) q[3];
sx q[3];
rz(2.5088333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.92155543) q[2];
sx q[2];
rz(-2.3735789) q[2];
sx q[2];
rz(1.7680232) q[2];
rz(0.31275648) q[3];
sx q[3];
rz(-2.0650568) q[3];
sx q[3];
rz(-3.055174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3304928) q[0];
sx q[0];
rz(-0.70527768) q[0];
sx q[0];
rz(0.16264859) q[0];
rz(1.9074408) q[1];
sx q[1];
rz(-1.994543) q[1];
sx q[1];
rz(2.2208234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4402078) q[0];
sx q[0];
rz(-0.21344859) q[0];
sx q[0];
rz(2.0646854) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92451325) q[2];
sx q[2];
rz(-1.1072031) q[2];
sx q[2];
rz(-1.4626423) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9879976) q[1];
sx q[1];
rz(-2.3210166) q[1];
sx q[1];
rz(1.6596036) q[1];
rz(0.7821857) q[3];
sx q[3];
rz(-1.4504004) q[3];
sx q[3];
rz(-0.52531717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2938701) q[2];
sx q[2];
rz(-0.72124481) q[2];
sx q[2];
rz(0.043962002) q[2];
rz(-1.7196767) q[3];
sx q[3];
rz(-0.43216643) q[3];
sx q[3];
rz(0.034984263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1304355) q[0];
sx q[0];
rz(-0.80428094) q[0];
sx q[0];
rz(-0.093611896) q[0];
rz(-1.0253819) q[1];
sx q[1];
rz(-1.5954115) q[1];
sx q[1];
rz(0.78790087) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2678515) q[0];
sx q[0];
rz(-2.2722167) q[0];
sx q[0];
rz(-2.3021477) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4008994) q[2];
sx q[2];
rz(-2.133103) q[2];
sx q[2];
rz(-3.135709) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.20847296) q[1];
sx q[1];
rz(-1.530024) q[1];
sx q[1];
rz(2.5426514) q[1];
rz(-2.2719194) q[3];
sx q[3];
rz(-0.18571412) q[3];
sx q[3];
rz(-2.7839422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0660144) q[2];
sx q[2];
rz(-0.54328537) q[2];
sx q[2];
rz(2.1755966) q[2];
rz(2.994359) q[3];
sx q[3];
rz(-1.4598264) q[3];
sx q[3];
rz(2.537263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0886993) q[0];
sx q[0];
rz(-1.3857144) q[0];
sx q[0];
rz(-1.8889486) q[0];
rz(0.057627536) q[1];
sx q[1];
rz(-1.3725955) q[1];
sx q[1];
rz(-2.7431814) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97038834) q[0];
sx q[0];
rz(-0.96231595) q[0];
sx q[0];
rz(1.7826016) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7407102) q[2];
sx q[2];
rz(-1.0771829) q[2];
sx q[2];
rz(2.5371671) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9908473) q[1];
sx q[1];
rz(-1.3708893) q[1];
sx q[1];
rz(-3.1374817) q[1];
rz(-pi) q[2];
rz(-2.5984405) q[3];
sx q[3];
rz(-1.4793605) q[3];
sx q[3];
rz(0.1426129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56732279) q[2];
sx q[2];
rz(-1.2719354) q[2];
sx q[2];
rz(2.2197913) q[2];
rz(0.70720339) q[3];
sx q[3];
rz(-2.0965529) q[3];
sx q[3];
rz(0.31064492) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1929753) q[0];
sx q[0];
rz(-3.0204168) q[0];
sx q[0];
rz(-2.2344672) q[0];
rz(-0.47997296) q[1];
sx q[1];
rz(-1.0555121) q[1];
sx q[1];
rz(-2.5352246) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9655749) q[0];
sx q[0];
rz(-1.2621573) q[0];
sx q[0];
rz(1.6762617) q[0];
x q[1];
rz(-2.1290413) q[2];
sx q[2];
rz(-1.1079567) q[2];
sx q[2];
rz(-1.3430581) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91057669) q[1];
sx q[1];
rz(-1.8635995) q[1];
sx q[1];
rz(-2.3412555) q[1];
rz(-pi) q[2];
rz(1.6799029) q[3];
sx q[3];
rz(-0.91489313) q[3];
sx q[3];
rz(0.65306585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2839462) q[2];
sx q[2];
rz(-1.5479156) q[2];
sx q[2];
rz(-2.471981) q[2];
rz(2.5745463) q[3];
sx q[3];
rz(-1.0457057) q[3];
sx q[3];
rz(-1.7414352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6794716) q[0];
sx q[0];
rz(-1.3422048) q[0];
sx q[0];
rz(0.8959499) q[0];
rz(0.040291928) q[1];
sx q[1];
rz(-2.1998684) q[1];
sx q[1];
rz(1.5641854) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5524254) q[0];
sx q[0];
rz(-1.5660813) q[0];
sx q[0];
rz(1.5943567) q[0];
rz(0.51367398) q[2];
sx q[2];
rz(-2.3875055) q[2];
sx q[2];
rz(0.23321345) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.40272199) q[1];
sx q[1];
rz(-0.67810692) q[1];
sx q[1];
rz(-0.88749591) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96503105) q[3];
sx q[3];
rz(-1.26767) q[3];
sx q[3];
rz(-2.3245963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3080052) q[2];
sx q[2];
rz(-0.31978017) q[2];
sx q[2];
rz(1.3099111) q[2];
rz(1.1854019) q[3];
sx q[3];
rz(-0.88806051) q[3];
sx q[3];
rz(1.6185224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0014521) q[0];
sx q[0];
rz(-1.3012713) q[0];
sx q[0];
rz(-0.95300994) q[0];
rz(0.078527191) q[1];
sx q[1];
rz(-1.086906) q[1];
sx q[1];
rz(-0.79798098) q[1];
rz(-2.2453515) q[2];
sx q[2];
rz(-1.7755388) q[2];
sx q[2];
rz(-0.15646738) q[2];
rz(-2.4918741) q[3];
sx q[3];
rz(-2.8856554) q[3];
sx q[3];
rz(1.9323991) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
