OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.69230429) q[0];
sx q[0];
rz(-0.46625724) q[0];
sx q[0];
rz(1.847108) q[0];
rz(-0.26588384) q[1];
sx q[1];
rz(-2.9335913) q[1];
sx q[1];
rz(1.2535569) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5781097) q[0];
sx q[0];
rz(-0.54395478) q[0];
sx q[0];
rz(-0.22814546) q[0];
rz(2.0090583) q[2];
sx q[2];
rz(-1.1992559) q[2];
sx q[2];
rz(-1.8105992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8575864) q[1];
sx q[1];
rz(-0.9269956) q[1];
sx q[1];
rz(2.8394152) q[1];
rz(-pi) q[2];
x q[2];
rz(0.01913602) q[3];
sx q[3];
rz(-1.2465256) q[3];
sx q[3];
rz(2.0024042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0365389) q[2];
sx q[2];
rz(-1.7221071) q[2];
sx q[2];
rz(-1.4744021) q[2];
rz(2.8640532) q[3];
sx q[3];
rz(-1.8614635) q[3];
sx q[3];
rz(-2.2138219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9934746) q[0];
sx q[0];
rz(-0.19522218) q[0];
sx q[0];
rz(1.8147722) q[0];
rz(-0.38816342) q[1];
sx q[1];
rz(-1.5048985) q[1];
sx q[1];
rz(0.51663748) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119902) q[0];
sx q[0];
rz(-1.4210295) q[0];
sx q[0];
rz(2.9277855) q[0];
rz(-pi) q[1];
rz(1.98968) q[2];
sx q[2];
rz(-1.3852714) q[2];
sx q[2];
rz(2.1806661) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.68143594) q[1];
sx q[1];
rz(-2.3189622) q[1];
sx q[1];
rz(2.9334738) q[1];
rz(-pi) q[2];
rz(-0.24833749) q[3];
sx q[3];
rz(-1.2988699) q[3];
sx q[3];
rz(0.14735315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8119729) q[2];
sx q[2];
rz(-0.5054349) q[2];
sx q[2];
rz(0.86414117) q[2];
rz(2.4498074) q[3];
sx q[3];
rz(-1.1312048) q[3];
sx q[3];
rz(2.2085021) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.5048051) q[1];
sx q[1];
rz(1.2791971) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7544781) q[0];
sx q[0];
rz(-0.25969782) q[0];
sx q[0];
rz(1.5952871) q[0];
x q[1];
rz(-2.6635936) q[2];
sx q[2];
rz(-2.1076116) q[2];
sx q[2];
rz(-1.9271242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0128263) q[1];
sx q[1];
rz(-0.44555659) q[1];
sx q[1];
rz(1.6548272) q[1];
rz(-0.62554977) q[3];
sx q[3];
rz(-1.8846785) q[3];
sx q[3];
rz(2.5108199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0731571) q[2];
sx q[2];
rz(-0.34978875) q[2];
sx q[2];
rz(-1.5901828) q[2];
rz(2.6326211) q[3];
sx q[3];
rz(-2.3059228) q[3];
sx q[3];
rz(-1.4510179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1125672) q[0];
sx q[0];
rz(-1.644716) q[0];
sx q[0];
rz(2.577884) q[0];
rz(-1.4504704) q[1];
sx q[1];
rz(-1.9861954) q[1];
sx q[1];
rz(2.3203826) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9552069) q[0];
sx q[0];
rz(-1.5816551) q[0];
sx q[0];
rz(-0.71594724) q[0];
x q[1];
rz(2.3963753) q[2];
sx q[2];
rz(-2.7181667) q[2];
sx q[2];
rz(-0.40027789) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5933696) q[1];
sx q[1];
rz(-2.0851622) q[1];
sx q[1];
rz(0.84131188) q[1];
x q[2];
rz(-0.037795732) q[3];
sx q[3];
rz(-2.1135847) q[3];
sx q[3];
rz(2.7181527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.700909) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(1.3779426) q[2];
rz(-2.9772229) q[3];
sx q[3];
rz(-1.5458958) q[3];
sx q[3];
rz(2.7730798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.565777) q[0];
sx q[0];
rz(-1.6933279) q[0];
sx q[0];
rz(-2.6337295) q[0];
rz(-0.85743088) q[1];
sx q[1];
rz(-1.6897759) q[1];
sx q[1];
rz(0.47971496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83905674) q[0];
sx q[0];
rz(-2.4094562) q[0];
sx q[0];
rz(-2.6468572) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7577997) q[2];
sx q[2];
rz(-1.6365882) q[2];
sx q[2];
rz(0.045595615) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.39998301) q[1];
sx q[1];
rz(-2.1711087) q[1];
sx q[1];
rz(0.27315477) q[1];
rz(-1.671319) q[3];
sx q[3];
rz(-2.8607603) q[3];
sx q[3];
rz(-2.7749244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.63626426) q[2];
sx q[2];
rz(-2.0365066) q[2];
sx q[2];
rz(2.9218033) q[2];
rz(0.2291186) q[3];
sx q[3];
rz(-0.90355211) q[3];
sx q[3];
rz(2.1598099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55766469) q[0];
sx q[0];
rz(-0.54323498) q[0];
sx q[0];
rz(2.014121) q[0];
rz(-0.30578956) q[1];
sx q[1];
rz(-1.1916279) q[1];
sx q[1];
rz(0.19010273) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33316225) q[0];
sx q[0];
rz(-1.8920533) q[0];
sx q[0];
rz(2.0684469) q[0];
rz(1.185955) q[2];
sx q[2];
rz(-2.3267496) q[2];
sx q[2];
rz(-0.84581918) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.99689647) q[1];
sx q[1];
rz(-2.3544925) q[1];
sx q[1];
rz(1.9793649) q[1];
rz(-pi) q[2];
rz(-2.8523731) q[3];
sx q[3];
rz(-2.8790124) q[3];
sx q[3];
rz(0.56910959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.019913435) q[2];
sx q[2];
rz(-1.6807669) q[2];
sx q[2];
rz(-2.0410247) q[2];
rz(1.8371643) q[3];
sx q[3];
rz(-1.5521939) q[3];
sx q[3];
rz(1.3531551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8562427) q[0];
sx q[0];
rz(-2.7142363) q[0];
sx q[0];
rz(0.57619488) q[0];
rz(-2.088749) q[1];
sx q[1];
rz(-0.57056999) q[1];
sx q[1];
rz(-2.0821234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5602901) q[0];
sx q[0];
rz(-1.6972739) q[0];
sx q[0];
rz(0.69445388) q[0];
rz(-pi) q[1];
rz(-3.1131831) q[2];
sx q[2];
rz(-1.1724262) q[2];
sx q[2];
rz(-2.0223126) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4814623) q[1];
sx q[1];
rz(-1.3795329) q[1];
sx q[1];
rz(0.44967117) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.430577) q[3];
sx q[3];
rz(-1.064015) q[3];
sx q[3];
rz(1.8975951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2759555) q[2];
sx q[2];
rz(-0.58791462) q[2];
sx q[2];
rz(2.0666583) q[2];
rz(-2.2771207) q[3];
sx q[3];
rz(-1.266022) q[3];
sx q[3];
rz(1.5629432) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7445755) q[0];
sx q[0];
rz(-0.96857849) q[0];
sx q[0];
rz(2.6343935) q[0];
rz(1.9604669) q[1];
sx q[1];
rz(-1.487251) q[1];
sx q[1];
rz(-0.91086737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024689704) q[0];
sx q[0];
rz(-2.0729613) q[0];
sx q[0];
rz(-0.85310081) q[0];
x q[1];
rz(2.9380625) q[2];
sx q[2];
rz(-0.99908057) q[2];
sx q[2];
rz(-2.2019486) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4994608) q[1];
sx q[1];
rz(-1.5593464) q[1];
sx q[1];
rz(1.4169372) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9796764) q[3];
sx q[3];
rz(-1.7954553) q[3];
sx q[3];
rz(2.4049559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18428093) q[0];
sx q[0];
rz(-0.98783699) q[0];
sx q[0];
rz(1.3405569) q[0];
rz(2.6181009) q[1];
sx q[1];
rz(-1.5308056) q[1];
sx q[1];
rz(-1.2045822) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8891994) q[0];
sx q[0];
rz(-2.3601818) q[0];
sx q[0];
rz(-1.8535421) q[0];
x q[1];
rz(0.14894375) q[2];
sx q[2];
rz(-1.7686378) q[2];
sx q[2];
rz(-0.94562519) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87540001) q[1];
sx q[1];
rz(-1.0169694) q[1];
sx q[1];
rz(-1.758105) q[1];
rz(-pi) q[2];
rz(-1.0597348) q[3];
sx q[3];
rz(-1.9283224) q[3];
sx q[3];
rz(-0.75393576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1022819) q[2];
sx q[2];
rz(-2.1151586) q[2];
sx q[2];
rz(-0.2956051) q[2];
rz(0.11232703) q[3];
sx q[3];
rz(-1.5887518) q[3];
sx q[3];
rz(2.2789392) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91117793) q[0];
sx q[0];
rz(-2.6417612) q[0];
sx q[0];
rz(2.3590132) q[0];
rz(-2.8453907) q[1];
sx q[1];
rz(-1.9007416) q[1];
sx q[1];
rz(0.20015073) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2252879) q[0];
sx q[0];
rz(-0.80919832) q[0];
sx q[0];
rz(-1.3455708) q[0];
rz(-0.28130071) q[2];
sx q[2];
rz(-2.1185115) q[2];
sx q[2];
rz(-0.74354625) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9505585) q[1];
sx q[1];
rz(-2.28015) q[1];
sx q[1];
rz(-0.85644958) q[1];
rz(-pi) q[2];
rz(2.2027722) q[3];
sx q[3];
rz(-2.0858313) q[3];
sx q[3];
rz(2.7523628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.062181648) q[2];
sx q[2];
rz(-0.14180413) q[2];
sx q[2];
rz(1.0395435) q[2];
rz(-3.1145575) q[3];
sx q[3];
rz(-2.1578433) q[3];
sx q[3];
rz(2.1496617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1886002) q[0];
sx q[0];
rz(-0.66467265) q[0];
sx q[0];
rz(-1.959214) q[0];
rz(1.4324808) q[1];
sx q[1];
rz(-1.2624337) q[1];
sx q[1];
rz(1.591325) q[1];
rz(-0.15122945) q[2];
sx q[2];
rz(-1.4561903) q[2];
sx q[2];
rz(-1.0753808) q[2];
rz(2.9441574) q[3];
sx q[3];
rz(-0.45396572) q[3];
sx q[3];
rz(2.288447) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
