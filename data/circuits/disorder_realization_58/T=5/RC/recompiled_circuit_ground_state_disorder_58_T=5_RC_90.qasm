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
rz(-1.2944846) q[0];
rz(-0.26588384) q[1];
sx q[1];
rz(-2.9335913) q[1];
sx q[1];
rz(1.2535569) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81119117) q[0];
sx q[0];
rz(-1.6881144) q[0];
sx q[0];
rz(0.53240029) q[0];
rz(0.40635477) q[2];
sx q[2];
rz(-1.1642805) q[2];
sx q[2];
rz(-2.7332566) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2840062) q[1];
sx q[1];
rz(-2.2145971) q[1];
sx q[1];
rz(2.8394152) q[1];
rz(0.01913602) q[3];
sx q[3];
rz(-1.2465256) q[3];
sx q[3];
rz(2.0024042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1050538) q[2];
sx q[2];
rz(-1.7221071) q[2];
sx q[2];
rz(1.6671906) q[2];
rz(2.8640532) q[3];
sx q[3];
rz(-1.2801291) q[3];
sx q[3];
rz(2.2138219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.14811806) q[0];
sx q[0];
rz(-0.19522218) q[0];
sx q[0];
rz(-1.3268205) q[0];
rz(-0.38816342) q[1];
sx q[1];
rz(-1.5048985) q[1];
sx q[1];
rz(0.51663748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119902) q[0];
sx q[0];
rz(-1.7205632) q[0];
sx q[0];
rz(-2.9277855) q[0];
rz(-pi) q[1];
rz(2.0031163) q[2];
sx q[2];
rz(-2.6856963) q[2];
sx q[2];
rz(2.1389463) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0320466) q[1];
sx q[1];
rz(-1.7228207) q[1];
sx q[1];
rz(0.8117453) q[1];
rz(-pi) q[2];
rz(-1.2907007) q[3];
sx q[3];
rz(-1.3317654) q[3];
sx q[3];
rz(1.7861507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3296198) q[2];
sx q[2];
rz(-2.6361578) q[2];
sx q[2];
rz(0.86414117) q[2];
rz(2.4498074) q[3];
sx q[3];
rz(-2.0103879) q[3];
sx q[3];
rz(-2.2085021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8253887) q[0];
sx q[0];
rz(-0.80556691) q[0];
sx q[0];
rz(-1.4027931) q[0];
rz(-2.5750776) q[1];
sx q[1];
rz(-1.6367876) q[1];
sx q[1];
rz(-1.8623955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3871146) q[0];
sx q[0];
rz(-0.25969782) q[0];
sx q[0];
rz(-1.5952871) q[0];
x q[1];
rz(-2.1612619) q[2];
sx q[2];
rz(-1.9771909) q[2];
sx q[2];
rz(-0.097336285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1058987) q[1];
sx q[1];
rz(-2.0146684) q[1];
sx q[1];
rz(0.040063362) q[1];
x q[2];
rz(-0.62554977) q[3];
sx q[3];
rz(-1.2569142) q[3];
sx q[3];
rz(-2.5108199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0731571) q[2];
sx q[2];
rz(-0.34978875) q[2];
sx q[2];
rz(-1.5901828) q[2];
rz(0.50897151) q[3];
sx q[3];
rz(-2.3059228) q[3];
sx q[3];
rz(-1.6905748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1125672) q[0];
sx q[0];
rz(-1.644716) q[0];
sx q[0];
rz(-2.577884) q[0];
rz(-1.6911223) q[1];
sx q[1];
rz(-1.9861954) q[1];
sx q[1];
rz(-2.3203826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7477363) q[0];
sx q[0];
rz(-0.85490037) q[0];
sx q[0];
rz(1.5564043) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8217373) q[2];
sx q[2];
rz(-1.8531688) q[2];
sx q[2];
rz(-2.6704463) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60821086) q[1];
sx q[1];
rz(-2.1898263) q[1];
sx q[1];
rz(2.4929895) q[1];
rz(-pi) q[2];
rz(1.508237) q[3];
sx q[3];
rz(-2.5976215) q[3];
sx q[3];
rz(0.49651745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.700909) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(-1.7636501) q[2];
rz(2.9772229) q[3];
sx q[3];
rz(-1.5458958) q[3];
sx q[3];
rz(-2.7730798) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57581562) q[0];
sx q[0];
rz(-1.6933279) q[0];
sx q[0];
rz(-2.6337295) q[0];
rz(-2.2841618) q[1];
sx q[1];
rz(-1.4518167) q[1];
sx q[1];
rz(0.47971496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83905674) q[0];
sx q[0];
rz(-2.4094562) q[0];
sx q[0];
rz(0.49473543) q[0];
rz(-1.4998593) q[2];
sx q[2];
rz(-1.1878769) q[2];
sx q[2];
rz(1.6429344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8138201) q[1];
sx q[1];
rz(-1.7952807) q[1];
sx q[1];
rz(-0.95275615) q[1];
rz(0.028939441) q[3];
sx q[3];
rz(-1.8501728) q[3];
sx q[3];
rz(-2.8795163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.63626426) q[2];
sx q[2];
rz(-1.1050861) q[2];
sx q[2];
rz(-2.9218033) q[2];
rz(0.2291186) q[3];
sx q[3];
rz(-0.90355211) q[3];
sx q[3];
rz(-0.98178274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.583928) q[0];
sx q[0];
rz(-0.54323498) q[0];
sx q[0];
rz(2.014121) q[0];
rz(-0.30578956) q[1];
sx q[1];
rz(-1.1916279) q[1];
sx q[1];
rz(-2.9514899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8084304) q[0];
sx q[0];
rz(-1.2495394) q[0];
sx q[0];
rz(-2.0684469) q[0];
rz(-1.185955) q[2];
sx q[2];
rz(-0.81484303) q[2];
sx q[2];
rz(-0.84581918) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2711004) q[1];
sx q[1];
rz(-1.8560579) q[1];
sx q[1];
rz(2.3149957) q[1];
rz(-pi) q[2];
x q[2];
rz(1.494287) q[3];
sx q[3];
rz(-1.8222295) q[3];
sx q[3];
rz(0.86800324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.019913435) q[2];
sx q[2];
rz(-1.4608258) q[2];
sx q[2];
rz(-2.0410247) q[2];
rz(-1.3044283) q[3];
sx q[3];
rz(-1.5521939) q[3];
sx q[3];
rz(-1.7884375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28534999) q[0];
sx q[0];
rz(-2.7142363) q[0];
sx q[0];
rz(-2.5653978) q[0];
rz(-1.0528437) q[1];
sx q[1];
rz(-0.57056999) q[1];
sx q[1];
rz(2.0821234) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.094166286) q[0];
sx q[0];
rz(-0.88297668) q[0];
sx q[0];
rz(1.4068013) q[0];
rz(-pi) q[1];
rz(1.6381864) q[2];
sx q[2];
rz(-2.7422649) q[2];
sx q[2];
rz(1.0461548) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53545836) q[1];
sx q[1];
rz(-2.6555057) q[1];
sx q[1];
rz(0.41907678) q[1];
rz(-2.430577) q[3];
sx q[3];
rz(-1.064015) q[3];
sx q[3];
rz(1.8975951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2759555) q[2];
sx q[2];
rz(-2.553678) q[2];
sx q[2];
rz(-1.0749344) q[2];
rz(-2.2771207) q[3];
sx q[3];
rz(-1.266022) q[3];
sx q[3];
rz(-1.5786494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39701715) q[0];
sx q[0];
rz(-2.1730142) q[0];
sx q[0];
rz(-2.6343935) q[0];
rz(-1.9604669) q[1];
sx q[1];
rz(-1.6543417) q[1];
sx q[1];
rz(2.2307253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0916101) q[0];
sx q[0];
rz(-2.292041) q[0];
sx q[0];
rz(0.87509416) q[0];
rz(-2.1520432) q[2];
sx q[2];
rz(-1.3999827) q[2];
sx q[2];
rz(0.51994158) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0720328) q[1];
sx q[1];
rz(-1.4169473) q[1];
sx q[1];
rz(-3.1300058) q[1];
rz(-pi) q[2];
rz(1.1619162) q[3];
sx q[3];
rz(-1.7954553) q[3];
sx q[3];
rz(-0.73663676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7562423) q[2];
sx q[2];
rz(-1.664398) q[2];
sx q[2];
rz(2.493646) q[2];
rz(-0.09662763) q[3];
sx q[3];
rz(-2.1164618) q[3];
sx q[3];
rz(2.1881762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18428093) q[0];
sx q[0];
rz(-2.1537557) q[0];
sx q[0];
rz(1.3405569) q[0];
rz(2.6181009) q[1];
sx q[1];
rz(-1.610787) q[1];
sx q[1];
rz(-1.9370105) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0055375) q[0];
sx q[0];
rz(-0.82804543) q[0];
sx q[0];
rz(0.27001794) q[0];
rz(-pi) q[1];
rz(2.2080293) q[2];
sx q[2];
rz(-2.8945403) q[2];
sx q[2];
rz(-2.8480094) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5455509) q[1];
sx q[1];
rz(-1.4117472) q[1];
sx q[1];
rz(-0.56175128) q[1];
rz(-pi) q[2];
rz(-0.91851652) q[3];
sx q[3];
rz(-2.5271086) q[3];
sx q[3];
rz(-2.8826734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.039310731) q[2];
sx q[2];
rz(-2.1151586) q[2];
sx q[2];
rz(0.2956051) q[2];
rz(-3.0292656) q[3];
sx q[3];
rz(-1.5528409) q[3];
sx q[3];
rz(0.86265341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91117793) q[0];
sx q[0];
rz(-2.6417612) q[0];
sx q[0];
rz(0.78257948) q[0];
rz(2.8453907) q[1];
sx q[1];
rz(-1.240851) q[1];
sx q[1];
rz(0.20015073) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3302932) q[0];
sx q[0];
rz(-1.7331373) q[0];
sx q[0];
rz(-0.77438023) q[0];
rz(-pi) q[1];
rz(-1.0050943) q[2];
sx q[2];
rz(-1.8100693) q[2];
sx q[2];
rz(0.67789652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0072277) q[1];
sx q[1];
rz(-1.0505465) q[1];
sx q[1];
rz(-2.2925329) q[1];
x q[2];
rz(-2.3347992) q[3];
sx q[3];
rz(-0.79232464) q[3];
sx q[3];
rz(-2.5522422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.062181648) q[2];
sx q[2];
rz(-2.9997885) q[2];
sx q[2];
rz(-2.1020491) q[2];
rz(-3.1145575) q[3];
sx q[3];
rz(-2.1578433) q[3];
sx q[3];
rz(2.1496617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95299245) q[0];
sx q[0];
rz(-2.47692) q[0];
sx q[0];
rz(1.1823786) q[0];
rz(-1.4324808) q[1];
sx q[1];
rz(-1.8791589) q[1];
sx q[1];
rz(-1.5502677) q[1];
rz(-2.9903632) q[2];
sx q[2];
rz(-1.6854023) q[2];
sx q[2];
rz(2.0662119) q[2];
rz(-1.6662206) q[3];
sx q[3];
rz(-2.0153042) q[3];
sx q[3];
rz(2.5074625) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
