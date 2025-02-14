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
rz(-2.6753354) q[0];
sx q[0];
rz(1.2944846) q[0];
rz(2.8757088) q[1];
sx q[1];
rz(2.9335913) q[1];
sx q[1];
rz(7.5367422) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3131375) q[0];
sx q[0];
rz(-2.0991517) q[0];
sx q[0];
rz(-1.4348475) q[0];
rz(-pi) q[1];
rz(-2.0090583) q[2];
sx q[2];
rz(-1.9423368) q[2];
sx q[2];
rz(-1.8105992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3365797) q[1];
sx q[1];
rz(-0.70195552) q[1];
sx q[1];
rz(-1.1932659) q[1];
rz(-pi) q[2];
rz(-3.1224566) q[3];
sx q[3];
rz(-1.8950671) q[3];
sx q[3];
rz(-2.0024042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1050538) q[2];
sx q[2];
rz(-1.4194856) q[2];
sx q[2];
rz(1.4744021) q[2];
rz(-2.8640532) q[3];
sx q[3];
rz(-1.2801291) q[3];
sx q[3];
rz(0.92777073) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14811806) q[0];
sx q[0];
rz(-2.9463705) q[0];
sx q[0];
rz(1.3268205) q[0];
rz(-0.38816342) q[1];
sx q[1];
rz(-1.5048985) q[1];
sx q[1];
rz(0.51663748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46088567) q[0];
sx q[0];
rz(-0.26038489) q[0];
sx q[0];
rz(0.61818122) q[0];
x q[1];
rz(-1.98968) q[2];
sx q[2];
rz(-1.3852714) q[2];
sx q[2];
rz(0.9609266) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.68143594) q[1];
sx q[1];
rz(-2.3189622) q[1];
sx q[1];
rz(-0.20811889) q[1];
x q[2];
rz(-2.8932552) q[3];
sx q[3];
rz(-1.8427227) q[3];
sx q[3];
rz(0.14735315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3296198) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8253887) q[0];
sx q[0];
rz(-0.80556691) q[0];
sx q[0];
rz(-1.7387996) q[0];
rz(2.5750776) q[1];
sx q[1];
rz(-1.6367876) q[1];
sx q[1];
rz(-1.2791971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7544781) q[0];
sx q[0];
rz(-2.8818948) q[0];
sx q[0];
rz(1.5952871) q[0];
rz(-pi) q[1];
rz(-2.1612619) q[2];
sx q[2];
rz(-1.9771909) q[2];
sx q[2];
rz(-0.097336285) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.5178901) q[1];
sx q[1];
rz(-1.5346171) q[1];
sx q[1];
rz(2.0149798) q[1];
x q[2];
rz(-0.62554977) q[3];
sx q[3];
rz(-1.2569142) q[3];
sx q[3];
rz(0.63077273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0684356) q[2];
sx q[2];
rz(-0.34978875) q[2];
sx q[2];
rz(1.5514099) q[2];
rz(0.50897151) q[3];
sx q[3];
rz(-2.3059228) q[3];
sx q[3];
rz(-1.6905748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.029025404) q[0];
sx q[0];
rz(-1.644716) q[0];
sx q[0];
rz(0.56370869) q[0];
rz(1.4504704) q[1];
sx q[1];
rz(-1.1553973) q[1];
sx q[1];
rz(-0.82121003) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37192908) q[0];
sx q[0];
rz(-0.71601501) q[0];
sx q[0];
rz(-0.016543702) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8674064) q[2];
sx q[2];
rz(-1.2640306) q[2];
sx q[2];
rz(-1.1916812) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5482231) q[1];
sx q[1];
rz(-2.0851622) q[1];
sx q[1];
rz(-2.3002808) q[1];
x q[2];
rz(-0.037795732) q[3];
sx q[3];
rz(-1.0280079) q[3];
sx q[3];
rz(-2.7181527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.700909) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(1.7636501) q[2];
rz(-0.16436973) q[3];
sx q[3];
rz(-1.5956968) q[3];
sx q[3];
rz(-0.36851287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.6618777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3025359) q[0];
sx q[0];
rz(-0.73213644) q[0];
sx q[0];
rz(2.6468572) q[0];
rz(-pi) q[1];
x q[1];
rz(0.38379294) q[2];
sx q[2];
rz(-1.6365882) q[2];
sx q[2];
rz(-3.095997) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8138201) q[1];
sx q[1];
rz(-1.7952807) q[1];
sx q[1];
rz(-2.1888365) q[1];
rz(-pi) q[2];
rz(-1.2913088) q[3];
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
rz(2.9218033) q[2];
rz(2.9124741) q[3];
sx q[3];
rz(-2.2380405) q[3];
sx q[3];
rz(-0.98178274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.583928) q[0];
sx q[0];
rz(-2.5983577) q[0];
sx q[0];
rz(-2.014121) q[0];
rz(0.30578956) q[1];
sx q[1];
rz(-1.9499648) q[1];
sx q[1];
rz(-2.9514899) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0738433) q[0];
sx q[0];
rz(-1.1007358) q[0];
sx q[0];
rz(2.7795634) q[0];
rz(-pi) q[1];
rz(-1.185955) q[2];
sx q[2];
rz(-0.81484303) q[2];
sx q[2];
rz(-0.84581918) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2711004) q[1];
sx q[1];
rz(-1.8560579) q[1];
sx q[1];
rz(-0.82659699) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8523731) q[3];
sx q[3];
rz(-2.8790124) q[3];
sx q[3];
rz(2.5724831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.019913435) q[2];
sx q[2];
rz(-1.6807669) q[2];
sx q[2];
rz(-2.0410247) q[2];
rz(-1.8371643) q[3];
sx q[3];
rz(-1.5893987) q[3];
sx q[3];
rz(-1.7884375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28534999) q[0];
sx q[0];
rz(-2.7142363) q[0];
sx q[0];
rz(0.57619488) q[0];
rz(-2.088749) q[1];
sx q[1];
rz(-2.5710227) q[1];
sx q[1];
rz(-1.0594692) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9807726) q[0];
sx q[0];
rz(-2.4376025) q[0];
sx q[0];
rz(2.9454548) q[0];
rz(-pi) q[1];
rz(-1.6381864) q[2];
sx q[2];
rz(-2.7422649) q[2];
sx q[2];
rz(-1.0461548) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66013038) q[1];
sx q[1];
rz(-1.3795329) q[1];
sx q[1];
rz(0.44967117) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93848159) q[3];
sx q[3];
rz(-2.1780067) q[3];
sx q[3];
rz(3.0724276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2759555) q[2];
sx q[2];
rz(-0.58791462) q[2];
sx q[2];
rz(-2.0666583) q[2];
rz(-2.2771207) q[3];
sx q[3];
rz(-1.8755707) q[3];
sx q[3];
rz(1.5786494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7445755) q[0];
sx q[0];
rz(-0.96857849) q[0];
sx q[0];
rz(0.50719914) q[0];
rz(1.9604669) q[1];
sx q[1];
rz(-1.6543417) q[1];
sx q[1];
rz(-2.2307253) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024689704) q[0];
sx q[0];
rz(-2.0729613) q[0];
sx q[0];
rz(2.2884918) q[0];
rz(2.9380625) q[2];
sx q[2];
rz(-2.1425121) q[2];
sx q[2];
rz(-0.93964404) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4994608) q[1];
sx q[1];
rz(-1.5593464) q[1];
sx q[1];
rz(-1.4169372) q[1];
x q[2];
rz(-2.0924545) q[3];
sx q[3];
rz(-0.46346617) q[3];
sx q[3];
rz(-1.3090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38535038) q[2];
sx q[2];
rz(-1.4771947) q[2];
sx q[2];
rz(-2.493646) q[2];
rz(0.09662763) q[3];
sx q[3];
rz(-1.0251309) q[3];
sx q[3];
rz(2.1881762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9573117) q[0];
sx q[0];
rz(-0.98783699) q[0];
sx q[0];
rz(-1.3405569) q[0];
rz(0.52349177) q[1];
sx q[1];
rz(-1.610787) q[1];
sx q[1];
rz(-1.2045822) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6197889) q[0];
sx q[0];
rz(-1.7685732) q[0];
sx q[0];
rz(2.3319582) q[0];
x q[1];
rz(-2.2080293) q[2];
sx q[2];
rz(-2.8945403) q[2];
sx q[2];
rz(-0.29358324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5455509) q[1];
sx q[1];
rz(-1.4117472) q[1];
sx q[1];
rz(-2.5798414) q[1];
rz(1.0597348) q[3];
sx q[3];
rz(-1.9283224) q[3];
sx q[3];
rz(-2.3876569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1022819) q[2];
sx q[2];
rz(-1.0264341) q[2];
sx q[2];
rz(2.8459876) q[2];
rz(-0.11232703) q[3];
sx q[3];
rz(-1.5528409) q[3];
sx q[3];
rz(2.2789392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91117793) q[0];
sx q[0];
rz(-0.4998315) q[0];
sx q[0];
rz(0.78257948) q[0];
rz(-0.29620194) q[1];
sx q[1];
rz(-1.240851) q[1];
sx q[1];
rz(0.20015073) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3302932) q[0];
sx q[0];
rz(-1.4084554) q[0];
sx q[0];
rz(-2.3672124) q[0];
x q[1];
rz(2.1364983) q[2];
sx q[2];
rz(-1.8100693) q[2];
sx q[2];
rz(0.67789652) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1175864) q[1];
sx q[1];
rz(-0.96029687) q[1];
sx q[1];
rz(2.4897051) q[1];
rz(2.3347992) q[3];
sx q[3];
rz(-0.79232464) q[3];
sx q[3];
rz(-0.58935048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.062181648) q[2];
sx q[2];
rz(-2.9997885) q[2];
sx q[2];
rz(-2.1020491) q[2];
rz(-0.027035106) q[3];
sx q[3];
rz(-2.1578433) q[3];
sx q[3];
rz(0.99193096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
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
rz(0.65244723) q[2];
sx q[2];
rz(-2.9521078) q[2];
sx q[2];
rz(-0.14833974) q[2];
rz(0.44628061) q[3];
sx q[3];
rz(-1.6569232) q[3];
sx q[3];
rz(0.89553064) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
