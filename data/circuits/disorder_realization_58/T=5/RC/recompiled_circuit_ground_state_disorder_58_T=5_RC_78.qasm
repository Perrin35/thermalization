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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5781097) q[0];
sx q[0];
rz(-0.54395478) q[0];
sx q[0];
rz(-0.22814546) q[0];
x q[1];
rz(-0.40635477) q[2];
sx q[2];
rz(-1.9773121) q[2];
sx q[2];
rz(0.40833607) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8575864) q[1];
sx q[1];
rz(-0.9269956) q[1];
sx q[1];
rz(0.30217742) q[1];
rz(-pi) q[2];
rz(1.8951224) q[3];
sx q[3];
rz(-1.5889349) q[3];
sx q[3];
rz(0.42551009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0365389) q[2];
sx q[2];
rz(-1.7221071) q[2];
sx q[2];
rz(-1.4744021) q[2];
rz(2.8640532) q[3];
sx q[3];
rz(-1.2801291) q[3];
sx q[3];
rz(-0.92777073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14811806) q[0];
sx q[0];
rz(-2.9463705) q[0];
sx q[0];
rz(-1.8147722) q[0];
rz(0.38816342) q[1];
sx q[1];
rz(-1.6366942) q[1];
sx q[1];
rz(-2.6249552) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.680707) q[0];
sx q[0];
rz(-0.26038489) q[0];
sx q[0];
rz(-2.5234114) q[0];
rz(-pi) q[1];
rz(1.98968) q[2];
sx q[2];
rz(-1.3852714) q[2];
sx q[2];
rz(2.1806661) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4601567) q[1];
sx q[1];
rz(-0.8226305) q[1];
sx q[1];
rz(0.20811889) q[1];
rz(-pi) q[2];
rz(-2.8932552) q[3];
sx q[3];
rz(-1.8427227) q[3];
sx q[3];
rz(-2.9942395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3296198) q[2];
sx q[2];
rz(-2.6361578) q[2];
sx q[2];
rz(2.2774515) q[2];
rz(-0.69178528) q[3];
sx q[3];
rz(-2.0103879) q[3];
sx q[3];
rz(0.9330906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8253887) q[0];
sx q[0];
rz(-0.80556691) q[0];
sx q[0];
rz(1.4027931) q[0];
rz(2.5750776) q[1];
sx q[1];
rz(-1.5048051) q[1];
sx q[1];
rz(1.2791971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3871146) q[0];
sx q[0];
rz(-0.25969782) q[0];
sx q[0];
rz(-1.5463056) q[0];
x q[1];
rz(-2.6635936) q[2];
sx q[2];
rz(-1.0339811) q[2];
sx q[2];
rz(-1.2144685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6237026) q[1];
sx q[1];
rz(-1.5346171) q[1];
sx q[1];
rz(-2.0149798) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1899105) q[3];
sx q[3];
rz(-0.98010585) q[3];
sx q[3];
rz(2.4209972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0731571) q[2];
sx q[2];
rz(-0.34978875) q[2];
sx q[2];
rz(1.5901828) q[2];
rz(-2.6326211) q[3];
sx q[3];
rz(-0.83566982) q[3];
sx q[3];
rz(-1.4510179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.029025404) q[0];
sx q[0];
rz(-1.4968766) q[0];
sx q[0];
rz(0.56370869) q[0];
rz(-1.6911223) q[1];
sx q[1];
rz(-1.9861954) q[1];
sx q[1];
rz(-2.3203826) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7696636) q[0];
sx q[0];
rz(-2.4255776) q[0];
sx q[0];
rz(-0.016543702) q[0];
x q[1];
rz(1.8674064) q[2];
sx q[2];
rz(-1.2640306) q[2];
sx q[2];
rz(1.1916812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.60821086) q[1];
sx q[1];
rz(-2.1898263) q[1];
sx q[1];
rz(0.64860313) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1037969) q[3];
sx q[3];
rz(-1.0280079) q[3];
sx q[3];
rz(0.42343994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.700909) q[2];
sx q[2];
rz(-2.4606885) q[2];
sx q[2];
rz(-1.7636501) q[2];
rz(-0.16436973) q[3];
sx q[3];
rz(-1.5956968) q[3];
sx q[3];
rz(-0.36851287) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57581562) q[0];
sx q[0];
rz(-1.4482647) q[0];
sx q[0];
rz(2.6337295) q[0];
rz(-0.85743088) q[1];
sx q[1];
rz(-1.6897759) q[1];
sx q[1];
rz(0.47971496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0282766) q[0];
sx q[0];
rz(-1.893766) q[0];
sx q[0];
rz(2.4723609) q[0];
rz(-1.6417333) q[2];
sx q[2];
rz(-1.9537158) q[2];
sx q[2];
rz(1.6429344) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7416096) q[1];
sx q[1];
rz(-0.97048391) q[1];
sx q[1];
rz(-2.8684379) q[1];
rz(-pi) q[2];
rz(1.2913088) q[3];
sx q[3];
rz(-1.5986134) q[3];
sx q[3];
rz(1.840855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5053284) q[2];
sx q[2];
rz(-1.1050861) q[2];
sx q[2];
rz(-2.9218033) q[2];
rz(-2.9124741) q[3];
sx q[3];
rz(-2.2380405) q[3];
sx q[3];
rz(0.98178274) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.583928) q[0];
sx q[0];
rz(-0.54323498) q[0];
sx q[0];
rz(1.1274717) q[0];
rz(-2.8358031) q[1];
sx q[1];
rz(-1.1916279) q[1];
sx q[1];
rz(2.9514899) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0738433) q[0];
sx q[0];
rz(-2.0408568) q[0];
sx q[0];
rz(-0.36202927) q[0];
rz(-pi) q[1];
rz(-0.79391329) q[2];
sx q[2];
rz(-1.2941235) q[2];
sx q[2];
rz(0.45396462) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87049228) q[1];
sx q[1];
rz(-1.8560579) q[1];
sx q[1];
rz(-2.3149957) q[1];
rz(-2.8894526) q[3];
sx q[3];
rz(-1.4966972) q[3];
sx q[3];
rz(-2.4197297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.019913435) q[2];
sx q[2];
rz(-1.4608258) q[2];
sx q[2];
rz(1.100568) q[2];
rz(-1.8371643) q[3];
sx q[3];
rz(-1.5893987) q[3];
sx q[3];
rz(1.3531551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28534999) q[0];
sx q[0];
rz(-0.42735639) q[0];
sx q[0];
rz(0.57619488) q[0];
rz(1.0528437) q[1];
sx q[1];
rz(-0.57056999) q[1];
sx q[1];
rz(-2.0821234) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9807726) q[0];
sx q[0];
rz(-2.4376025) q[0];
sx q[0];
rz(-2.9454548) q[0];
x q[1];
rz(0.028409516) q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-0.66013038) q[1];
sx q[1];
rz(-1.3795329) q[1];
sx q[1];
rz(2.6919215) q[1];
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
x q[1];
rz(1.8656371) q[2];
sx q[2];
rz(-2.553678) q[2];
sx q[2];
rz(-2.0666583) q[2];
rz(-2.2771207) q[3];
sx q[3];
rz(-1.266022) q[3];
sx q[3];
rz(1.5629432) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7445755) q[0];
sx q[0];
rz(-2.1730142) q[0];
sx q[0];
rz(-2.6343935) q[0];
rz(-1.9604669) q[1];
sx q[1];
rz(-1.487251) q[1];
sx q[1];
rz(-2.2307253) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024689704) q[0];
sx q[0];
rz(-2.0729613) q[0];
sx q[0];
rz(-0.85310081) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2664001) q[2];
sx q[2];
rz(-2.5385529) q[2];
sx q[2];
rz(-1.3040486) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14503059) q[1];
sx q[1];
rz(-2.9873114) q[1];
sx q[1];
rz(1.645374) q[1];
x q[2];
rz(0.24408079) q[3];
sx q[3];
rz(-1.968813) q[3];
sx q[3];
rz(0.73792968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38535038) q[2];
sx q[2];
rz(-1.4771947) q[2];
sx q[2];
rz(2.493646) q[2];
rz(3.044965) q[3];
sx q[3];
rz(-1.0251309) q[3];
sx q[3];
rz(0.9534165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.18428093) q[0];
sx q[0];
rz(-0.98783699) q[0];
sx q[0];
rz(1.3405569) q[0];
rz(-0.52349177) q[1];
sx q[1];
rz(-1.5308056) q[1];
sx q[1];
rz(1.9370105) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25239326) q[0];
sx q[0];
rz(-0.78141087) q[0];
sx q[0];
rz(1.8535421) q[0];
rz(-pi) q[1];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-0.87540001) q[1];
sx q[1];
rz(-2.1246233) q[1];
sx q[1];
rz(1.3834877) q[1];
rz(-pi) q[2];
rz(-0.91851652) q[3];
sx q[3];
rz(-2.5271086) q[3];
sx q[3];
rz(0.25891925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1022819) q[2];
sx q[2];
rz(-2.1151586) q[2];
sx q[2];
rz(-2.8459876) q[2];
rz(-3.0292656) q[3];
sx q[3];
rz(-1.5528409) q[3];
sx q[3];
rz(-2.2789392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91117793) q[0];
sx q[0];
rz(-2.6417612) q[0];
sx q[0];
rz(-0.78257948) q[0];
rz(-0.29620194) q[1];
sx q[1];
rz(-1.9007416) q[1];
sx q[1];
rz(-0.20015073) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59574612) q[0];
sx q[0];
rz(-0.78773088) q[0];
sx q[0];
rz(-2.911522) q[0];
x q[1];
rz(0.28130071) q[2];
sx q[2];
rz(-2.1185115) q[2];
sx q[2];
rz(0.74354625) q[2];
rz(-pi) q[3];
x q[3];
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
rz(-0.65188758) q[1];
rz(-pi) q[2];
rz(-2.2027722) q[3];
sx q[3];
rz(-1.0557613) q[3];
sx q[3];
rz(-0.38922986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.079411) q[2];
sx q[2];
rz(-2.9997885) q[2];
sx q[2];
rz(-2.1020491) q[2];
rz(0.027035106) q[3];
sx q[3];
rz(-2.1578433) q[3];
sx q[3];
rz(-0.99193096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.1886002) q[0];
sx q[0];
rz(-2.47692) q[0];
sx q[0];
rz(1.1823786) q[0];
rz(1.7091119) q[1];
sx q[1];
rz(-1.8791589) q[1];
sx q[1];
rz(-1.5502677) q[1];
rz(-1.4548789) q[2];
sx q[2];
rz(-1.4205665) q[2];
sx q[2];
rz(0.51284075) q[2];
rz(2.695312) q[3];
sx q[3];
rz(-1.4846694) q[3];
sx q[3];
rz(-2.246062) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
