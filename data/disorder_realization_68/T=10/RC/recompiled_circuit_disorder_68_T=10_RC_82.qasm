OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7473937) q[0];
sx q[0];
rz(-2.6497901) q[0];
sx q[0];
rz(2.9536182) q[0];
rz(2.0239053) q[1];
sx q[1];
rz(4.6586577) q[1];
sx q[1];
rz(12.933856) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1023094) q[0];
sx q[0];
rz(-1.4674205) q[0];
sx q[0];
rz(-2.1084059) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0796702) q[2];
sx q[2];
rz(-1.6899781) q[2];
sx q[2];
rz(-2.116938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.623466) q[1];
sx q[1];
rz(-1.7382442) q[1];
sx q[1];
rz(-1.3256339) q[1];
rz(-pi) q[2];
rz(0.34100451) q[3];
sx q[3];
rz(-1.738428) q[3];
sx q[3];
rz(1.0238907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1774896) q[2];
sx q[2];
rz(-2.6314645) q[2];
sx q[2];
rz(-2.5906079) q[2];
rz(1.3059113) q[3];
sx q[3];
rz(-1.6492313) q[3];
sx q[3];
rz(-1.3163542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6630163) q[0];
sx q[0];
rz(-1.0000279) q[0];
sx q[0];
rz(2.6696894) q[0];
rz(-0.42981237) q[1];
sx q[1];
rz(-1.8919573) q[1];
sx q[1];
rz(-0.93634161) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72915709) q[0];
sx q[0];
rz(-1.658174) q[0];
sx q[0];
rz(-0.23075128) q[0];
rz(-pi) q[1];
rz(-0.41994862) q[2];
sx q[2];
rz(-2.311085) q[2];
sx q[2];
rz(-0.74479693) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4221103) q[1];
sx q[1];
rz(-1.6750095) q[1];
sx q[1];
rz(2.4476493) q[1];
rz(-pi) q[2];
rz(1.068088) q[3];
sx q[3];
rz(-1.5503251) q[3];
sx q[3];
rz(-2.3308144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.77461809) q[2];
sx q[2];
rz(-2.8145511) q[2];
sx q[2];
rz(0.42638391) q[2];
rz(-1.9042227) q[3];
sx q[3];
rz(-2.5137413) q[3];
sx q[3];
rz(0.0330851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8957829) q[0];
sx q[0];
rz(-1.824546) q[0];
sx q[0];
rz(-2.202503) q[0];
rz(-2.242873) q[1];
sx q[1];
rz(-2.6627314) q[1];
sx q[1];
rz(-0.59392196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9868601) q[0];
sx q[0];
rz(-1.8911456) q[0];
sx q[0];
rz(0.31517584) q[0];
x q[1];
rz(-2.2592696) q[2];
sx q[2];
rz(-1.5261298) q[2];
sx q[2];
rz(-2.4633212) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9582639) q[1];
sx q[1];
rz(-1.3576344) q[1];
sx q[1];
rz(1.1413241) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46843723) q[3];
sx q[3];
rz(-1.2611946) q[3];
sx q[3];
rz(0.34164159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.64017355) q[2];
sx q[2];
rz(-0.72945014) q[2];
sx q[2];
rz(-1.4397941) q[2];
rz(2.7539608) q[3];
sx q[3];
rz(-1.6250316) q[3];
sx q[3];
rz(0.38813996) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3751635) q[0];
sx q[0];
rz(-1.5901934) q[0];
sx q[0];
rz(-2.638812) q[0];
rz(-0.76820961) q[1];
sx q[1];
rz(-2.6380824) q[1];
sx q[1];
rz(-2.3847413) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4810281) q[0];
sx q[0];
rz(-1.9303983) q[0];
sx q[0];
rz(-1.2147551) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5966162) q[2];
sx q[2];
rz(-0.46586793) q[2];
sx q[2];
rz(-2.4398838) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4591603) q[1];
sx q[1];
rz(-1.2130514) q[1];
sx q[1];
rz(2.1898502) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5960856) q[3];
sx q[3];
rz(-2.3295998) q[3];
sx q[3];
rz(0.10520392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7148774) q[2];
sx q[2];
rz(-1.9116414) q[2];
sx q[2];
rz(1.654401) q[2];
rz(2.5590844) q[3];
sx q[3];
rz(-2.0472066) q[3];
sx q[3];
rz(0.55707651) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476167) q[0];
sx q[0];
rz(-2.0563545) q[0];
sx q[0];
rz(0.75772444) q[0];
rz(1.853653) q[1];
sx q[1];
rz(-0.92823354) q[1];
sx q[1];
rz(1.0505189) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5245363) q[0];
sx q[0];
rz(-0.40818383) q[0];
sx q[0];
rz(0.88390669) q[0];
rz(-2.1804785) q[2];
sx q[2];
rz(-2.4498307) q[2];
sx q[2];
rz(-1.320653) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.98002906) q[1];
sx q[1];
rz(-1.0161576) q[1];
sx q[1];
rz(2.5142923) q[1];
rz(-pi) q[2];
rz(1.0195144) q[3];
sx q[3];
rz(-0.84421221) q[3];
sx q[3];
rz(-2.7725078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.22333764) q[2];
sx q[2];
rz(-0.35327062) q[2];
sx q[2];
rz(0.63344947) q[2];
rz(1.194362) q[3];
sx q[3];
rz(-1.4712237) q[3];
sx q[3];
rz(-2.4244394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.441992) q[0];
sx q[0];
rz(-2.6514335) q[0];
sx q[0];
rz(2.8884086) q[0];
rz(1.5340012) q[1];
sx q[1];
rz(-1.4350767) q[1];
sx q[1];
rz(1.6794499) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77308649) q[0];
sx q[0];
rz(-2.9162772) q[0];
sx q[0];
rz(-0.36264514) q[0];
rz(-0.24121933) q[2];
sx q[2];
rz(-2.2058645) q[2];
sx q[2];
rz(-1.6068174) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4436795) q[1];
sx q[1];
rz(-1.4005449) q[1];
sx q[1];
rz(2.7531569) q[1];
rz(-2.8263894) q[3];
sx q[3];
rz(-1.8990371) q[3];
sx q[3];
rz(2.4344276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.2016466) q[2];
sx q[2];
rz(-2.3942409) q[2];
sx q[2];
rz(-0.80491006) q[2];
rz(-1.1770052) q[3];
sx q[3];
rz(-2.1046808) q[3];
sx q[3];
rz(-1.0837519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8686304) q[0];
sx q[0];
rz(-1.0694163) q[0];
sx q[0];
rz(-2.2139363) q[0];
rz(-2.1169128) q[1];
sx q[1];
rz(-1.506348) q[1];
sx q[1];
rz(-2.129508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9818078) q[0];
sx q[0];
rz(-0.73043981) q[0];
sx q[0];
rz(-2.3083789) q[0];
rz(-pi) q[1];
rz(0.38883932) q[2];
sx q[2];
rz(-1.6992339) q[2];
sx q[2];
rz(-1.625979) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5092897) q[1];
sx q[1];
rz(-1.0273233) q[1];
sx q[1];
rz(-1.5322881) q[1];
rz(-pi) q[2];
rz(-2.5829671) q[3];
sx q[3];
rz(-1.6405676) q[3];
sx q[3];
rz(0.39486265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6107789) q[2];
sx q[2];
rz(-1.4799708) q[2];
sx q[2];
rz(-0.42993316) q[2];
rz(-1.0144462) q[3];
sx q[3];
rz(-0.40922624) q[3];
sx q[3];
rz(0.53340069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6124509) q[0];
sx q[0];
rz(-0.93944678) q[0];
sx q[0];
rz(2.9274143) q[0];
rz(2.0902436) q[1];
sx q[1];
rz(-2.9290757) q[1];
sx q[1];
rz(0.28373757) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8003214) q[0];
sx q[0];
rz(-2.8386136) q[0];
sx q[0];
rz(-0.11462258) q[0];
rz(-pi) q[1];
rz(2.0869414) q[2];
sx q[2];
rz(-0.37545855) q[2];
sx q[2];
rz(1.6106538) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.46854308) q[1];
sx q[1];
rz(-2.0646411) q[1];
sx q[1];
rz(-1.7935351) q[1];
rz(0.65225668) q[3];
sx q[3];
rz(-1.0271003) q[3];
sx q[3];
rz(0.039747681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6909137) q[2];
sx q[2];
rz(-0.51270715) q[2];
sx q[2];
rz(-1.696375) q[2];
rz(-1.5444267) q[3];
sx q[3];
rz(-1.4404567) q[3];
sx q[3];
rz(2.8022695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63672367) q[0];
sx q[0];
rz(-0.050014194) q[0];
sx q[0];
rz(-0.069256393) q[0];
rz(-1.4878558) q[1];
sx q[1];
rz(-1.8585049) q[1];
sx q[1];
rz(1.5725296) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5208961) q[0];
sx q[0];
rz(-0.81917742) q[0];
sx q[0];
rz(-2.2197414) q[0];
rz(0.32585085) q[2];
sx q[2];
rz(-1.648765) q[2];
sx q[2];
rz(-0.98380145) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.442765) q[1];
sx q[1];
rz(-2.8932533) q[1];
sx q[1];
rz(0.34766867) q[1];
rz(-pi) q[2];
rz(0.2089573) q[3];
sx q[3];
rz(-2.2088802) q[3];
sx q[3];
rz(-1.0554505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9562324) q[2];
sx q[2];
rz(-2.9113443) q[2];
sx q[2];
rz(3.0017079) q[2];
rz(-0.36758962) q[3];
sx q[3];
rz(-1.9544173) q[3];
sx q[3];
rz(0.99115133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1763828) q[0];
sx q[0];
rz(-0.39127025) q[0];
sx q[0];
rz(2.4998253) q[0];
rz(1.9104674) q[1];
sx q[1];
rz(-1.9893913) q[1];
sx q[1];
rz(-0.26783255) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3788911) q[0];
sx q[0];
rz(-0.91103103) q[0];
sx q[0];
rz(-0.99878175) q[0];
x q[1];
rz(-1.7540625) q[2];
sx q[2];
rz(-1.6943309) q[2];
sx q[2];
rz(2.0545517) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1416727) q[1];
sx q[1];
rz(-1.281764) q[1];
sx q[1];
rz(-2.8029867) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.855741) q[3];
sx q[3];
rz(-1.0470069) q[3];
sx q[3];
rz(-1.7459735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.81007593) q[2];
sx q[2];
rz(-1.6167567) q[2];
sx q[2];
rz(-1.4769185) q[2];
rz(-0.26633513) q[3];
sx q[3];
rz(-0.24644066) q[3];
sx q[3];
rz(2.5951071) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.01263604) q[0];
sx q[0];
rz(-0.92100443) q[0];
sx q[0];
rz(-2.2367649) q[0];
rz(2.3616882) q[1];
sx q[1];
rz(-0.48702469) q[1];
sx q[1];
rz(-1.3866966) q[1];
rz(-0.98942479) q[2];
sx q[2];
rz(-2.1944254) q[2];
sx q[2];
rz(2.2133322) q[2];
rz(-2.4640502) q[3];
sx q[3];
rz(-1.4031706) q[3];
sx q[3];
rz(-1.8085898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
