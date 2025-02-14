OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4492884) q[0];
sx q[0];
rz(3.6078499) q[0];
sx q[0];
rz(10.719263) q[0];
rz(2.8757088) q[1];
sx q[1];
rz(-0.20800132) q[1];
sx q[1];
rz(1.8880358) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3131375) q[0];
sx q[0];
rz(-2.0991517) q[0];
sx q[0];
rz(-1.7067451) q[0];
rz(-pi) q[1];
rz(0.40635477) q[2];
sx q[2];
rz(-1.9773121) q[2];
sx q[2];
rz(-0.40833607) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2840062) q[1];
sx q[1];
rz(-0.9269956) q[1];
sx q[1];
rz(0.30217742) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5139318) q[3];
sx q[3];
rz(-0.32481501) q[3];
sx q[3];
rz(-1.9424094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0365389) q[2];
sx q[2];
rz(-1.4194856) q[2];
sx q[2];
rz(1.6671906) q[2];
rz(-2.8640532) q[3];
sx q[3];
rz(-1.2801291) q[3];
sx q[3];
rz(-2.2138219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14811806) q[0];
sx q[0];
rz(-0.19522218) q[0];
sx q[0];
rz(1.3268205) q[0];
rz(-2.7534292) q[1];
sx q[1];
rz(-1.6366942) q[1];
sx q[1];
rz(0.51663748) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1735794) q[0];
sx q[0];
rz(-1.359419) q[0];
sx q[0];
rz(-1.7239991) q[0];
rz(-pi) q[1];
rz(-0.20262589) q[2];
sx q[2];
rz(-1.1595402) q[2];
sx q[2];
rz(-2.613668) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7611446) q[1];
sx q[1];
rz(-0.7711322) q[1];
sx q[1];
rz(-1.7898331) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2932986) q[3];
sx q[3];
rz(-2.775421) q[3];
sx q[3];
rz(2.2375905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3296198) q[2];
sx q[2];
rz(-2.6361578) q[2];
sx q[2];
rz(2.2774515) q[2];
rz(-2.4498074) q[3];
sx q[3];
rz(-2.0103879) q[3];
sx q[3];
rz(-0.9330906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.316204) q[0];
sx q[0];
rz(-0.80556691) q[0];
sx q[0];
rz(-1.7387996) q[0];
rz(2.5750776) q[1];
sx q[1];
rz(-1.6367876) q[1];
sx q[1];
rz(1.8623955) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7544781) q[0];
sx q[0];
rz(-2.8818948) q[0];
sx q[0];
rz(-1.5463056) q[0];
x q[1];
rz(0.91275349) q[2];
sx q[2];
rz(-0.70281723) q[2];
sx q[2];
rz(1.1352486) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0128263) q[1];
sx q[1];
rz(-2.6960361) q[1];
sx q[1];
rz(-1.6548272) q[1];
x q[2];
rz(0.50619979) q[3];
sx q[3];
rz(-2.4512614) q[3];
sx q[3];
rz(-1.3439646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0731571) q[2];
sx q[2];
rz(-2.7918039) q[2];
sx q[2];
rz(-1.5901828) q[2];
rz(0.50897151) q[3];
sx q[3];
rz(-2.3059228) q[3];
sx q[3];
rz(1.4510179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029025404) q[0];
sx q[0];
rz(-1.644716) q[0];
sx q[0];
rz(-2.577884) q[0];
rz(-1.6911223) q[1];
sx q[1];
rz(-1.9861954) q[1];
sx q[1];
rz(0.82121003) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1863857) q[0];
sx q[0];
rz(-1.5816551) q[0];
sx q[0];
rz(-0.71594724) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3963753) q[2];
sx q[2];
rz(-2.7181667) q[2];
sx q[2];
rz(2.7413148) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5933696) q[1];
sx q[1];
rz(-1.0564305) q[1];
sx q[1];
rz(-2.3002808) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0276919) q[3];
sx q[3];
rz(-1.6031577) q[3];
sx q[3];
rz(-2.0137656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.700909) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(1.3779426) q[2];
rz(0.16436973) q[3];
sx q[3];
rz(-1.5458958) q[3];
sx q[3];
rz(2.7730798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57581562) q[0];
sx q[0];
rz(-1.4482647) q[0];
sx q[0];
rz(-0.50786316) q[0];
rz(-0.85743088) q[1];
sx q[1];
rz(-1.4518167) q[1];
sx q[1];
rz(-0.47971496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21151152) q[0];
sx q[0];
rz(-2.1997612) q[0];
sx q[0];
rz(-1.9741365) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4998593) q[2];
sx q[2];
rz(-1.1878769) q[2];
sx q[2];
rz(1.6429344) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39998301) q[1];
sx q[1];
rz(-2.1711087) q[1];
sx q[1];
rz(2.8684379) q[1];
rz(-pi) q[2];
rz(1.8502838) q[3];
sx q[3];
rz(-1.5986134) q[3];
sx q[3];
rz(1.3007377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5053284) q[2];
sx q[2];
rz(-1.1050861) q[2];
sx q[2];
rz(-0.21978933) q[2];
rz(-0.2291186) q[3];
sx q[3];
rz(-0.90355211) q[3];
sx q[3];
rz(-2.1598099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.583928) q[0];
sx q[0];
rz(-0.54323498) q[0];
sx q[0];
rz(-2.014121) q[0];
rz(2.8358031) q[1];
sx q[1];
rz(-1.1916279) q[1];
sx q[1];
rz(-2.9514899) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33316225) q[0];
sx q[0];
rz(-1.8920533) q[0];
sx q[0];
rz(-2.0684469) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9556376) q[2];
sx q[2];
rz(-2.3267496) q[2];
sx q[2];
rz(-2.2957735) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87049228) q[1];
sx q[1];
rz(-1.2855347) q[1];
sx q[1];
rz(-2.3149957) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25214002) q[3];
sx q[3];
rz(-1.4966972) q[3];
sx q[3];
rz(0.72186294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.019913435) q[2];
sx q[2];
rz(-1.6807669) q[2];
sx q[2];
rz(-1.100568) q[2];
rz(-1.8371643) q[3];
sx q[3];
rz(-1.5893987) q[3];
sx q[3];
rz(1.3531551) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28534999) q[0];
sx q[0];
rz(-0.42735639) q[0];
sx q[0];
rz(-2.5653978) q[0];
rz(-1.0528437) q[1];
sx q[1];
rz(-2.5710227) q[1];
sx q[1];
rz(1.0594692) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094166286) q[0];
sx q[0];
rz(-0.88297668) q[0];
sx q[0];
rz(1.7347914) q[0];
rz(-pi) q[1];
rz(1.5034063) q[2];
sx q[2];
rz(-2.7422649) q[2];
sx q[2];
rz(-1.0461548) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66013038) q[1];
sx q[1];
rz(-1.3795329) q[1];
sx q[1];
rz(-0.44967117) q[1];
rz(-pi) q[2];
rz(-0.70487421) q[3];
sx q[3];
rz(-0.84669152) q[3];
sx q[3];
rz(2.3016229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8656371) q[2];
sx q[2];
rz(-0.58791462) q[2];
sx q[2];
rz(-2.0666583) q[2];
rz(2.2771207) q[3];
sx q[3];
rz(-1.266022) q[3];
sx q[3];
rz(1.5786494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39701715) q[0];
sx q[0];
rz(-0.96857849) q[0];
sx q[0];
rz(-2.6343935) q[0];
rz(1.1811258) q[1];
sx q[1];
rz(-1.6543417) q[1];
sx q[1];
rz(-0.91086737) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0499826) q[0];
sx q[0];
rz(-2.292041) q[0];
sx q[0];
rz(0.87509416) q[0];
rz(-pi) q[1];
rz(-2.1520432) q[2];
sx q[2];
rz(-1.74161) q[2];
sx q[2];
rz(-0.51994158) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.14503059) q[1];
sx q[1];
rz(-0.15428126) q[1];
sx q[1];
rz(1.645374) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0924545) q[3];
sx q[3];
rz(-0.46346617) q[3];
sx q[3];
rz(1.3090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7562423) q[2];
sx q[2];
rz(-1.4771947) q[2];
sx q[2];
rz(0.64794668) q[2];
rz(3.044965) q[3];
sx q[3];
rz(-1.0251309) q[3];
sx q[3];
rz(-2.1881762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9573117) q[0];
sx q[0];
rz(-0.98783699) q[0];
sx q[0];
rz(-1.8010358) q[0];
rz(-2.6181009) q[1];
sx q[1];
rz(-1.610787) q[1];
sx q[1];
rz(1.9370105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6197889) q[0];
sx q[0];
rz(-1.3730195) q[0];
sx q[0];
rz(-0.80963442) q[0];
rz(-pi) q[1];
rz(-2.9926489) q[2];
sx q[2];
rz(-1.7686378) q[2];
sx q[2];
rz(2.1959675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9203176) q[1];
sx q[1];
rz(-2.560096) q[1];
sx q[1];
rz(0.29249549) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0597348) q[3];
sx q[3];
rz(-1.2132702) q[3];
sx q[3];
rz(2.3876569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1022819) q[2];
sx q[2];
rz(-2.1151586) q[2];
sx q[2];
rz(-2.8459876) q[2];
rz(-3.0292656) q[3];
sx q[3];
rz(-1.5887518) q[3];
sx q[3];
rz(-0.86265341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2304147) q[0];
sx q[0];
rz(-0.4998315) q[0];
sx q[0];
rz(0.78257948) q[0];
rz(-0.29620194) q[1];
sx q[1];
rz(-1.240851) q[1];
sx q[1];
rz(0.20015073) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91630477) q[0];
sx q[0];
rz(-2.3323943) q[0];
sx q[0];
rz(-1.3455708) q[0];
rz(1.0050943) q[2];
sx q[2];
rz(-1.8100693) q[2];
sx q[2];
rz(-0.67789652) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.134365) q[1];
sx q[1];
rz(-1.0505465) q[1];
sx q[1];
rz(2.2925329) q[1];
rz(-pi) q[2];
rz(0.61171054) q[3];
sx q[3];
rz(-2.1107622) q[3];
sx q[3];
rz(-1.5276791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.079411) q[2];
sx q[2];
rz(-2.9997885) q[2];
sx q[2];
rz(-1.0395435) q[2];
rz(-0.027035106) q[3];
sx q[3];
rz(-2.1578433) q[3];
sx q[3];
rz(0.99193096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1886002) q[0];
sx q[0];
rz(-2.47692) q[0];
sx q[0];
rz(1.1823786) q[0];
rz(-1.7091119) q[1];
sx q[1];
rz(-1.2624337) q[1];
sx q[1];
rz(1.591325) q[1];
rz(1.6867137) q[2];
sx q[2];
rz(-1.4205665) q[2];
sx q[2];
rz(0.51284075) q[2];
rz(-1.4753721) q[3];
sx q[3];
rz(-1.1262885) q[3];
sx q[3];
rz(-0.63413017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
