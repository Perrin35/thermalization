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
rz(0.82062757) q[0];
sx q[0];
rz(-0.66250116) q[0];
sx q[0];
rz(2.878046) q[0];
rz(-2.0343556) q[1];
sx q[1];
rz(-1.2821481) q[1];
sx q[1];
rz(2.4660304) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5476108) q[0];
sx q[0];
rz(-1.2351961) q[0];
sx q[0];
rz(-0.006860126) q[0];
rz(-pi) q[1];
rz(2.1525429) q[2];
sx q[2];
rz(-2.0768696) q[2];
sx q[2];
rz(-2.842234) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1141211) q[1];
sx q[1];
rz(-0.41931376) q[1];
sx q[1];
rz(1.9472576) q[1];
rz(-pi) q[2];
x q[2];
rz(0.3278424) q[3];
sx q[3];
rz(-2.050542) q[3];
sx q[3];
rz(0.2827417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.88966173) q[2];
sx q[2];
rz(-1.9193005) q[2];
sx q[2];
rz(0.3678073) q[2];
rz(0.31428549) q[3];
sx q[3];
rz(-0.29044423) q[3];
sx q[3];
rz(0.17816003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80938584) q[0];
sx q[0];
rz(-3.0520913) q[0];
sx q[0];
rz(1.7820763) q[0];
rz(1.8571732) q[1];
sx q[1];
rz(-1.9764683) q[1];
sx q[1];
rz(-0.051609106) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4415539) q[0];
sx q[0];
rz(-1.7518028) q[0];
sx q[0];
rz(-1.3010398) q[0];
x q[1];
rz(-1.2003635) q[2];
sx q[2];
rz(-1.7631754) q[2];
sx q[2];
rz(2.7967601) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7588531) q[1];
sx q[1];
rz(-1.5713646) q[1];
sx q[1];
rz(1.9447536) q[1];
rz(-0.69689023) q[3];
sx q[3];
rz(-2.5514388) q[3];
sx q[3];
rz(-1.3499638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.08903683) q[2];
sx q[2];
rz(-2.8812228) q[2];
sx q[2];
rz(-0.54711771) q[2];
rz(-1.268092) q[3];
sx q[3];
rz(-1.3601466) q[3];
sx q[3];
rz(-0.33856302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7608305) q[0];
sx q[0];
rz(-2.1050504) q[0];
sx q[0];
rz(1.3408252) q[0];
rz(-2.2876168) q[1];
sx q[1];
rz(-1.8546591) q[1];
sx q[1];
rz(0.30240789) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9895565) q[0];
sx q[0];
rz(-1.8468509) q[0];
sx q[0];
rz(-0.92854519) q[0];
x q[1];
rz(0.45291169) q[2];
sx q[2];
rz(-2.6476417) q[2];
sx q[2];
rz(-0.043722186) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21061347) q[1];
sx q[1];
rz(-0.96081464) q[1];
sx q[1];
rz(0.90175866) q[1];
x q[2];
rz(2.0539927) q[3];
sx q[3];
rz(-0.54387605) q[3];
sx q[3];
rz(0.15009201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2053908) q[2];
sx q[2];
rz(-1.256029) q[2];
sx q[2];
rz(2.5993247) q[2];
rz(-0.99644709) q[3];
sx q[3];
rz(-0.22765972) q[3];
sx q[3];
rz(-3.0136133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
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
rz(0.19313136) q[0];
sx q[0];
rz(-1.6574991) q[0];
sx q[0];
rz(1.8782072) q[0];
rz(-0.43426934) q[1];
sx q[1];
rz(-2.3691005) q[1];
sx q[1];
rz(-1.6901406) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49187351) q[0];
sx q[0];
rz(-1.8611835) q[0];
sx q[0];
rz(0.47652737) q[0];
rz(-0.66777467) q[2];
sx q[2];
rz(-0.66168565) q[2];
sx q[2];
rz(-2.0354605) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8122708) q[1];
sx q[1];
rz(-1.4854317) q[1];
sx q[1];
rz(-1.4567767) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2882388) q[3];
sx q[3];
rz(-2.382016) q[3];
sx q[3];
rz(-0.64379287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0666535) q[2];
sx q[2];
rz(-0.94330072) q[2];
sx q[2];
rz(-2.2554876) q[2];
rz(-3.1352299) q[3];
sx q[3];
rz(-1.3402091) q[3];
sx q[3];
rz(1.1191012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56115443) q[0];
sx q[0];
rz(-0.4564603) q[0];
sx q[0];
rz(1.2503257) q[0];
rz(-0.23288947) q[1];
sx q[1];
rz(-2.3857375) q[1];
sx q[1];
rz(-0.09045352) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0165774) q[0];
sx q[0];
rz(-2.8629486) q[0];
sx q[0];
rz(-2.0296627) q[0];
rz(-1.301342) q[2];
sx q[2];
rz(-1.1134992) q[2];
sx q[2];
rz(-2.2895165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4793746) q[1];
sx q[1];
rz(-1.5081128) q[1];
sx q[1];
rz(-2.4880243) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8188263) q[3];
sx q[3];
rz(-2.8042386) q[3];
sx q[3];
rz(0.017291822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6470486) q[2];
sx q[2];
rz(-0.19379751) q[2];
sx q[2];
rz(-1.2079283) q[2];
rz(0.9211933) q[3];
sx q[3];
rz(-2.2974206) q[3];
sx q[3];
rz(-2.6827961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0334466) q[0];
sx q[0];
rz(-1.7195846) q[0];
sx q[0];
rz(0.49149996) q[0];
rz(0.72832251) q[1];
sx q[1];
rz(-1.1110543) q[1];
sx q[1];
rz(-0.19457766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96724961) q[0];
sx q[0];
rz(-1.3493269) q[0];
sx q[0];
rz(-0.99507777) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75677432) q[2];
sx q[2];
rz(-2.1770453) q[2];
sx q[2];
rz(2.221125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.814333) q[1];
sx q[1];
rz(-2.2047298) q[1];
sx q[1];
rz(-2.4933706) q[1];
rz(-1.3380906) q[3];
sx q[3];
rz(-0.92706087) q[3];
sx q[3];
rz(0.88730592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7877833) q[2];
sx q[2];
rz(-2.1259191) q[2];
sx q[2];
rz(0.0018399012) q[2];
rz(1.4932102) q[3];
sx q[3];
rz(-2.4108678) q[3];
sx q[3];
rz(-2.8572237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2581185) q[0];
sx q[0];
rz(-2.879877) q[0];
sx q[0];
rz(2.0765685) q[0];
rz(0.26271543) q[1];
sx q[1];
rz(-2.5698667) q[1];
sx q[1];
rz(2.6782742) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7119448) q[0];
sx q[0];
rz(-1.4139626) q[0];
sx q[0];
rz(-2.1615087) q[0];
x q[1];
rz(-2.1049785) q[2];
sx q[2];
rz(-0.91832238) q[2];
sx q[2];
rz(-0.62803113) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83144298) q[1];
sx q[1];
rz(-2.0528704) q[1];
sx q[1];
rz(-1.7245913) q[1];
x q[2];
rz(3.1172905) q[3];
sx q[3];
rz(-1.6858471) q[3];
sx q[3];
rz(2.3772511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.00039438417) q[2];
sx q[2];
rz(-1.980915) q[2];
sx q[2];
rz(0.82994962) q[2];
rz(2.0937894) q[3];
sx q[3];
rz(-2.6959097) q[3];
sx q[3];
rz(-2.9629663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8038427) q[0];
sx q[0];
rz(-0.99030817) q[0];
sx q[0];
rz(-0.75482279) q[0];
rz(2.7524475) q[1];
sx q[1];
rz(-1.6837589) q[1];
sx q[1];
rz(0.39774242) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11618092) q[0];
sx q[0];
rz(-1.8278024) q[0];
sx q[0];
rz(-1.8639131) q[0];
x q[1];
rz(-0.36093485) q[2];
sx q[2];
rz(-0.58141469) q[2];
sx q[2];
rz(-1.6045398) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1017572) q[1];
sx q[1];
rz(-1.3164742) q[1];
sx q[1];
rz(-1.8524342) q[1];
x q[2];
rz(-2.4598148) q[3];
sx q[3];
rz(-2.2684134) q[3];
sx q[3];
rz(0.57428065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8336746) q[2];
sx q[2];
rz(-0.29611823) q[2];
sx q[2];
rz(0.31416848) q[2];
rz(1.8999772) q[3];
sx q[3];
rz(-1.5528468) q[3];
sx q[3];
rz(1.9023021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.817713) q[0];
sx q[0];
rz(-2.4535024) q[0];
sx q[0];
rz(0.24539465) q[0];
rz(-1.9112401) q[1];
sx q[1];
rz(-0.10086682) q[1];
sx q[1];
rz(-2.6255677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2986261) q[0];
sx q[0];
rz(-0.43382513) q[0];
sx q[0];
rz(-2.6515657) q[0];
rz(0.7528968) q[2];
sx q[2];
rz(-2.12924) q[2];
sx q[2];
rz(0.073592984) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5799595) q[1];
sx q[1];
rz(-1.206534) q[1];
sx q[1];
rz(1.2666525) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13652674) q[3];
sx q[3];
rz(-1.6322005) q[3];
sx q[3];
rz(-1.9207973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5838991) q[2];
sx q[2];
rz(-1.2217174) q[2];
sx q[2];
rz(-1.0830797) q[2];
rz(2.9205868) q[3];
sx q[3];
rz(-2.998304) q[3];
sx q[3];
rz(2.3451282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021127064) q[0];
sx q[0];
rz(-3.0607304) q[0];
sx q[0];
rz(3.0057111) q[0];
rz(1.7120301) q[1];
sx q[1];
rz(-2.0202961) q[1];
sx q[1];
rz(-0.28724614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82086876) q[0];
sx q[0];
rz(-2.831376) q[0];
sx q[0];
rz(-0.91281548) q[0];
rz(-0.084566074) q[2];
sx q[2];
rz(-1.8075602) q[2];
sx q[2];
rz(2.6344476) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.72037) q[1];
sx q[1];
rz(-2.1541641) q[1];
sx q[1];
rz(-1.3066584) q[1];
rz(-pi) q[2];
rz(1.997825) q[3];
sx q[3];
rz(-0.62814071) q[3];
sx q[3];
rz(-2.028156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2222432) q[2];
sx q[2];
rz(-1.9025849) q[2];
sx q[2];
rz(-0.83402056) q[2];
rz(-2.833994) q[3];
sx q[3];
rz(-0.37720507) q[3];
sx q[3];
rz(2.8789177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.96497) q[0];
sx q[0];
rz(-1.8360092) q[0];
sx q[0];
rz(-1.0259884) q[0];
rz(2.7550244) q[1];
sx q[1];
rz(-1.3810806) q[1];
sx q[1];
rz(2.5234533) q[1];
rz(-1.7025399) q[2];
sx q[2];
rz(-1.9921039) q[2];
sx q[2];
rz(2.8767246) q[2];
rz(-1.781842) q[3];
sx q[3];
rz(-1.7132218) q[3];
sx q[3];
rz(1.6794459) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
