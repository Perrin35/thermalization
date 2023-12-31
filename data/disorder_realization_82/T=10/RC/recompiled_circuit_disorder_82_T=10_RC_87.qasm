OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80469552) q[0];
sx q[0];
rz(-1.0372294) q[0];
sx q[0];
rz(0.35559911) q[0];
rz(-0.30272499) q[1];
sx q[1];
rz(-2.0974789) q[1];
sx q[1];
rz(1.8619327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049392603) q[0];
sx q[0];
rz(-2.85673) q[0];
sx q[0];
rz(-2.0869135) q[0];
x q[1];
rz(1.7945292) q[2];
sx q[2];
rz(-1.1587843) q[2];
sx q[2];
rz(-1.7681233) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0919839) q[1];
sx q[1];
rz(-2.1503452) q[1];
sx q[1];
rz(-1.0773354) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9079307) q[3];
sx q[3];
rz(-1.7710925) q[3];
sx q[3];
rz(1.3113126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1352284) q[2];
sx q[2];
rz(-0.93868119) q[2];
sx q[2];
rz(-0.65650666) q[2];
rz(-0.73900977) q[3];
sx q[3];
rz(-2.6813172) q[3];
sx q[3];
rz(0.41729331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0083369) q[0];
sx q[0];
rz(-2.3454741) q[0];
sx q[0];
rz(-2.546229) q[0];
rz(3.0796675) q[1];
sx q[1];
rz(-1.9134816) q[1];
sx q[1];
rz(-0.48746902) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66747626) q[0];
sx q[0];
rz(-1.5063018) q[0];
sx q[0];
rz(-1.6518031) q[0];
x q[1];
rz(-2.7198118) q[2];
sx q[2];
rz(-1.9739208) q[2];
sx q[2];
rz(0.031907206) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0812379) q[1];
sx q[1];
rz(-1.3831257) q[1];
sx q[1];
rz(-0.16138046) q[1];
x q[2];
rz(3.0549166) q[3];
sx q[3];
rz(-1.8868586) q[3];
sx q[3];
rz(-0.3294657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.1797103) q[2];
sx q[2];
rz(-1.8790745) q[2];
sx q[2];
rz(-2.7462192) q[2];
rz(1.0428492) q[3];
sx q[3];
rz(-0.53111774) q[3];
sx q[3];
rz(-0.038671967) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19668002) q[0];
sx q[0];
rz(-2.0356464) q[0];
sx q[0];
rz(2.9512067) q[0];
rz(3.0186675) q[1];
sx q[1];
rz(-0.38750896) q[1];
sx q[1];
rz(2.9188459) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4091464) q[0];
sx q[0];
rz(-1.1799066) q[0];
sx q[0];
rz(1.976165) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0388971) q[2];
sx q[2];
rz(-0.42429081) q[2];
sx q[2];
rz(-1.588856) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.45914868) q[1];
sx q[1];
rz(-0.61668452) q[1];
sx q[1];
rz(3.0241806) q[1];
rz(-pi) q[2];
rz(0.95007105) q[3];
sx q[3];
rz(-0.40682236) q[3];
sx q[3];
rz(1.0328968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8824076) q[2];
sx q[2];
rz(-1.8683445) q[2];
sx q[2];
rz(-1.9474691) q[2];
rz(-1.0549226) q[3];
sx q[3];
rz(-1.6566365) q[3];
sx q[3];
rz(2.1779493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38917437) q[0];
sx q[0];
rz(-0.27802935) q[0];
sx q[0];
rz(1.7383204) q[0];
rz(1.6482884) q[1];
sx q[1];
rz(-1.5416668) q[1];
sx q[1];
rz(0.9202252) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7677778) q[0];
sx q[0];
rz(-2.2245363) q[0];
sx q[0];
rz(-0.52315229) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.780519) q[2];
sx q[2];
rz(-1.4898584) q[2];
sx q[2];
rz(1.0635478) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8060382) q[1];
sx q[1];
rz(-2.028855) q[1];
sx q[1];
rz(2.4469417) q[1];
rz(2.782015) q[3];
sx q[3];
rz(-1.7633071) q[3];
sx q[3];
rz(-1.0958835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.197864) q[2];
sx q[2];
rz(-1.4131763) q[2];
sx q[2];
rz(1.8614004) q[2];
rz(0.81930339) q[3];
sx q[3];
rz(-1.3179444) q[3];
sx q[3];
rz(1.7839446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.458805) q[0];
sx q[0];
rz(-1.7651432) q[0];
sx q[0];
rz(1.4047594) q[0];
rz(2.3732896) q[1];
sx q[1];
rz(-0.65015018) q[1];
sx q[1];
rz(0.4531025) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4015409) q[0];
sx q[0];
rz(-1.6443559) q[0];
sx q[0];
rz(-3.1042276) q[0];
rz(-pi) q[1];
rz(-1.4278973) q[2];
sx q[2];
rz(-1.60534) q[2];
sx q[2];
rz(2.3847716) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.15935005) q[1];
sx q[1];
rz(-2.9061926) q[1];
sx q[1];
rz(1.9297286) q[1];
rz(0.46053967) q[3];
sx q[3];
rz(-2.1201049) q[3];
sx q[3];
rz(2.3159546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.12864628) q[2];
sx q[2];
rz(-2.3036028) q[2];
sx q[2];
rz(-1.5117234) q[2];
rz(-2.417918) q[3];
sx q[3];
rz(-1.2440163) q[3];
sx q[3];
rz(0.85030142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.033427514) q[0];
sx q[0];
rz(-1.7280248) q[0];
sx q[0];
rz(2.9034555) q[0];
rz(-2.761633) q[1];
sx q[1];
rz(-2.0894876) q[1];
sx q[1];
rz(-1.628081) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3753525) q[0];
sx q[0];
rz(-1.1002812) q[0];
sx q[0];
rz(1.9348295) q[0];
rz(-pi) q[1];
rz(1.6983301) q[2];
sx q[2];
rz(-2.6886534) q[2];
sx q[2];
rz(-2.4085101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5781128) q[1];
sx q[1];
rz(-1.1155177) q[1];
sx q[1];
rz(-0.77225765) q[1];
x q[2];
rz(2.2699039) q[3];
sx q[3];
rz(-2.5786434) q[3];
sx q[3];
rz(2.9339919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0825519) q[2];
sx q[2];
rz(-2.5566792) q[2];
sx q[2];
rz(1.9980105) q[2];
rz(-0.13051662) q[3];
sx q[3];
rz(-1.7139707) q[3];
sx q[3];
rz(0.17091621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0156353) q[0];
sx q[0];
rz(-0.97706777) q[0];
sx q[0];
rz(-0.39757279) q[0];
rz(-1.3145087) q[1];
sx q[1];
rz(-0.54388261) q[1];
sx q[1];
rz(-1.1869173) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2448954) q[0];
sx q[0];
rz(-2.5573686) q[0];
sx q[0];
rz(2.0969735) q[0];
rz(-pi) q[1];
x q[1];
rz(0.059139472) q[2];
sx q[2];
rz(-1.5416317) q[2];
sx q[2];
rz(-1.189122) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0901511) q[1];
sx q[1];
rz(-1.7045583) q[1];
sx q[1];
rz(-0.83530463) q[1];
rz(2.5116634) q[3];
sx q[3];
rz(-1.808872) q[3];
sx q[3];
rz(-1.5743953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.879803) q[2];
sx q[2];
rz(-1.3061085) q[2];
sx q[2];
rz(-1.3558033) q[2];
rz(1.5073744) q[3];
sx q[3];
rz(-0.58871388) q[3];
sx q[3];
rz(-0.9986977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81925201) q[0];
sx q[0];
rz(-1.9544019) q[0];
sx q[0];
rz(-0.85574714) q[0];
rz(-3.1198655) q[1];
sx q[1];
rz(-2.1179492) q[1];
sx q[1];
rz(-1.0303248) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54661575) q[0];
sx q[0];
rz(-2.3175276) q[0];
sx q[0];
rz(1.0975518) q[0];
x q[1];
rz(-0.048642283) q[2];
sx q[2];
rz(-2.0970793) q[2];
sx q[2];
rz(0.7359879) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9955666) q[1];
sx q[1];
rz(-1.4194173) q[1];
sx q[1];
rz(1.8712908) q[1];
x q[2];
rz(0.52328531) q[3];
sx q[3];
rz(-2.0998294) q[3];
sx q[3];
rz(-2.7494591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8490303) q[2];
sx q[2];
rz(-1.8886671) q[2];
sx q[2];
rz(-0.070177468) q[2];
rz(2.3146546) q[3];
sx q[3];
rz(-1.6707872) q[3];
sx q[3];
rz(-2.5551445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4940015) q[0];
sx q[0];
rz(-1.4275455) q[0];
sx q[0];
rz(2.8253187) q[0];
rz(2.0896185) q[1];
sx q[1];
rz(-0.72967044) q[1];
sx q[1];
rz(-2.0064328) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90056706) q[0];
sx q[0];
rz(-1.3824438) q[0];
sx q[0];
rz(-2.7450949) q[0];
rz(0.029149292) q[2];
sx q[2];
rz(-1.8910732) q[2];
sx q[2];
rz(-2.0707891) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.70442048) q[1];
sx q[1];
rz(-1.6295027) q[1];
sx q[1];
rz(1.0874332) q[1];
x q[2];
rz(1.6039861) q[3];
sx q[3];
rz(-0.6233218) q[3];
sx q[3];
rz(0.25191307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6241374) q[2];
sx q[2];
rz(-0.74666658) q[2];
sx q[2];
rz(1.5967782) q[2];
rz(0.67772135) q[3];
sx q[3];
rz(-0.8422519) q[3];
sx q[3];
rz(0.28963447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.233376) q[0];
sx q[0];
rz(-0.69013086) q[0];
sx q[0];
rz(2.7897575) q[0];
rz(-2.8219163) q[1];
sx q[1];
rz(-0.37477481) q[1];
sx q[1];
rz(0.19616729) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10712121) q[0];
sx q[0];
rz(-2.2959318) q[0];
sx q[0];
rz(0.79191533) q[0];
rz(-1.4354544) q[2];
sx q[2];
rz(-1.5208828) q[2];
sx q[2];
rz(2.4382255) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4539459) q[1];
sx q[1];
rz(-1.5554264) q[1];
sx q[1];
rz(-0.83666283) q[1];
rz(-pi) q[2];
rz(-0.84368002) q[3];
sx q[3];
rz(-1.0495249) q[3];
sx q[3];
rz(-0.47714864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3907884) q[2];
sx q[2];
rz(-1.5020341) q[2];
sx q[2];
rz(0.081136726) q[2];
rz(1.1817415) q[3];
sx q[3];
rz(-0.90098444) q[3];
sx q[3];
rz(2.0666163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(-0.022973013) q[0];
sx q[0];
rz(-1.58283) q[0];
sx q[0];
rz(1.3118623) q[0];
rz(-1.3148057) q[1];
sx q[1];
rz(-0.79827764) q[1];
sx q[1];
rz(-1.6006443) q[1];
rz(-2.9076004) q[2];
sx q[2];
rz(-1.6178314) q[2];
sx q[2];
rz(-2.9754054) q[2];
rz(3.0951981) q[3];
sx q[3];
rz(-1.3413324) q[3];
sx q[3];
rz(-3.0975292) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
