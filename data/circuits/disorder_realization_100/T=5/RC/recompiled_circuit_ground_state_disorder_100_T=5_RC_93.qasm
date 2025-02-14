OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.72803175) q[0];
sx q[0];
rz(2.1821238) q[0];
sx q[0];
rz(9.9766599) q[0];
rz(-0.38504398) q[1];
sx q[1];
rz(-1.3621962) q[1];
sx q[1];
rz(0.41866067) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075567632) q[0];
sx q[0];
rz(-2.2097144) q[0];
sx q[0];
rz(0.35982168) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7483226) q[2];
sx q[2];
rz(-1.1273618) q[2];
sx q[2];
rz(0.73278945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6578411) q[1];
sx q[1];
rz(-0.22210177) q[1];
sx q[1];
rz(-0.41390837) q[1];
rz(-1.6012264) q[3];
sx q[3];
rz(-1.55369) q[3];
sx q[3];
rz(0.78029437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5358413) q[2];
sx q[2];
rz(-1.302482) q[2];
sx q[2];
rz(-2.7525986) q[2];
rz(2.244106) q[3];
sx q[3];
rz(-2.5806081) q[3];
sx q[3];
rz(1.2787904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9101343) q[0];
sx q[0];
rz(-2.1899905) q[0];
sx q[0];
rz(-2.8125473) q[0];
rz(1.344205) q[1];
sx q[1];
rz(-0.75733328) q[1];
sx q[1];
rz(-0.12408852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0253925) q[0];
sx q[0];
rz(-1.3611462) q[0];
sx q[0];
rz(1.8127182) q[0];
x q[1];
rz(-0.28606881) q[2];
sx q[2];
rz(-1.7033615) q[2];
sx q[2];
rz(2.2372467) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8491148) q[1];
sx q[1];
rz(-0.58365959) q[1];
sx q[1];
rz(0.8846585) q[1];
rz(-0.20601087) q[3];
sx q[3];
rz(-2.2643746) q[3];
sx q[3];
rz(1.3242974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.2825534) q[2];
sx q[2];
rz(-2.6586847) q[2];
sx q[2];
rz(-2.5340951) q[2];
rz(2.2533158) q[3];
sx q[3];
rz(-1.6162623) q[3];
sx q[3];
rz(-1.7150778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68179503) q[0];
sx q[0];
rz(-1.8495704) q[0];
sx q[0];
rz(0.23455308) q[0];
rz(-2.081743) q[1];
sx q[1];
rz(-1.1666965) q[1];
sx q[1];
rz(-2.2134728) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9938894) q[0];
sx q[0];
rz(-0.77489955) q[0];
sx q[0];
rz(2.5323917) q[0];
rz(-pi) q[1];
rz(1.1546722) q[2];
sx q[2];
rz(-1.1319931) q[2];
sx q[2];
rz(-0.051953944) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4473572) q[1];
sx q[1];
rz(-2.7856845) q[1];
sx q[1];
rz(0.39519542) q[1];
rz(1.9506128) q[3];
sx q[3];
rz(-0.39576021) q[3];
sx q[3];
rz(2.2896965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5782535) q[2];
sx q[2];
rz(-1.4207062) q[2];
sx q[2];
rz(1.8948179) q[2];
rz(0.16683821) q[3];
sx q[3];
rz(-0.92883795) q[3];
sx q[3];
rz(1.5290574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.700915) q[0];
sx q[0];
rz(-0.066078521) q[0];
sx q[0];
rz(-2.3833158) q[0];
rz(-1.5658763) q[1];
sx q[1];
rz(-2.5570452) q[1];
sx q[1];
rz(2.27104) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1681002) q[0];
sx q[0];
rz(-2.3273558) q[0];
sx q[0];
rz(-0.6310985) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0921368) q[2];
sx q[2];
rz(-1.0564197) q[2];
sx q[2];
rz(-2.5829022) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4728001) q[1];
sx q[1];
rz(-1.973098) q[1];
sx q[1];
rz(-2.975906) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4985256) q[3];
sx q[3];
rz(-2.1974583) q[3];
sx q[3];
rz(-2.5570803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42252758) q[2];
sx q[2];
rz(-2.1336522) q[2];
sx q[2];
rz(-0.75971216) q[2];
rz(-0.64940137) q[3];
sx q[3];
rz(-1.5348744) q[3];
sx q[3];
rz(-1.5627129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0287057) q[0];
sx q[0];
rz(-0.6539456) q[0];
sx q[0];
rz(1.9586067) q[0];
rz(-0.68823367) q[1];
sx q[1];
rz(-1.3166683) q[1];
sx q[1];
rz(2.7722955) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0205905) q[0];
sx q[0];
rz(-0.90743104) q[0];
sx q[0];
rz(-1.536002) q[0];
x q[1];
rz(2.4299116) q[2];
sx q[2];
rz(-1.3698915) q[2];
sx q[2];
rz(2.7679659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5253882) q[1];
sx q[1];
rz(-1.3862351) q[1];
sx q[1];
rz(-2.2730458) q[1];
rz(0.33325382) q[3];
sx q[3];
rz(-1.8849533) q[3];
sx q[3];
rz(-0.21018782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1987622) q[2];
sx q[2];
rz(-1.8161215) q[2];
sx q[2];
rz(2.14373) q[2];
rz(2.5241847) q[3];
sx q[3];
rz(-2.182775) q[3];
sx q[3];
rz(-0.82120419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1029516) q[0];
sx q[0];
rz(-1.2741673) q[0];
sx q[0];
rz(-1.2432903) q[0];
rz(0.78168166) q[1];
sx q[1];
rz(-1.2739173) q[1];
sx q[1];
rz(-2.8807358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2170625) q[0];
sx q[0];
rz(-2.6017761) q[0];
sx q[0];
rz(0.41882078) q[0];
rz(-pi) q[1];
rz(0.37520295) q[2];
sx q[2];
rz(-1.979036) q[2];
sx q[2];
rz(0.21438504) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1697537) q[1];
sx q[1];
rz(-2.0287645) q[1];
sx q[1];
rz(2.9724246) q[1];
rz(-pi) q[2];
rz(-1.8571768) q[3];
sx q[3];
rz(-2.4900511) q[3];
sx q[3];
rz(-2.1809354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47238749) q[2];
sx q[2];
rz(-0.86296764) q[2];
sx q[2];
rz(-0.48259398) q[2];
rz(-0.11416301) q[3];
sx q[3];
rz(-1.0897021) q[3];
sx q[3];
rz(-2.193006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.362185) q[0];
sx q[0];
rz(-3.0857093) q[0];
sx q[0];
rz(0.26595297) q[0];
rz(-2.7159122) q[1];
sx q[1];
rz(-1.8961779) q[1];
sx q[1];
rz(0.46806213) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47129318) q[0];
sx q[0];
rz(-2.6001996) q[0];
sx q[0];
rz(-1.5575717) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8276026) q[2];
sx q[2];
rz(-2.7903284) q[2];
sx q[2];
rz(-2.6333269) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0056161) q[1];
sx q[1];
rz(-1.8274725) q[1];
sx q[1];
rz(-2.0784432) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1105177) q[3];
sx q[3];
rz(-1.0508176) q[3];
sx q[3];
rz(-1.7620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.257306) q[2];
sx q[2];
rz(-1.851119) q[2];
sx q[2];
rz(-2.5878944) q[2];
rz(0.92343679) q[3];
sx q[3];
rz(-2.2348576) q[3];
sx q[3];
rz(-3.0807909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4629352) q[0];
sx q[0];
rz(-0.23660062) q[0];
sx q[0];
rz(-0.67805725) q[0];
rz(-2.1642115) q[1];
sx q[1];
rz(-1.1481608) q[1];
sx q[1];
rz(-1.703702) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2861833) q[0];
sx q[0];
rz(-1.1146422) q[0];
sx q[0];
rz(2.9265704) q[0];
rz(-pi) q[1];
rz(-1.5735658) q[2];
sx q[2];
rz(-1.0068147) q[2];
sx q[2];
rz(-2.0071047) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9919093) q[1];
sx q[1];
rz(-1.339773) q[1];
sx q[1];
rz(2.8808589) q[1];
rz(0.74681654) q[3];
sx q[3];
rz(-1.1594611) q[3];
sx q[3];
rz(-2.8432027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4108654) q[2];
sx q[2];
rz(-2.5812456) q[2];
sx q[2];
rz(-2.4746258) q[2];
rz(-0.11735958) q[3];
sx q[3];
rz(-1.2753692) q[3];
sx q[3];
rz(-1.3807152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9475107) q[0];
sx q[0];
rz(-1.5503333) q[0];
sx q[0];
rz(-0.34341735) q[0];
rz(1.6471479) q[1];
sx q[1];
rz(-0.98973715) q[1];
sx q[1];
rz(2.6572878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5268516) q[0];
sx q[0];
rz(-2.1526045) q[0];
sx q[0];
rz(1.2528166) q[0];
rz(2.7126409) q[2];
sx q[2];
rz(-1.9520734) q[2];
sx q[2];
rz(-2.5457053) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.81496325) q[1];
sx q[1];
rz(-2.1479534) q[1];
sx q[1];
rz(1.3054791) q[1];
rz(1.7913489) q[3];
sx q[3];
rz(-1.327264) q[3];
sx q[3];
rz(1.1597553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9453498) q[2];
sx q[2];
rz(-2.0413155) q[2];
sx q[2];
rz(2.2958882) q[2];
rz(-1.9218933) q[3];
sx q[3];
rz(-1.2518576) q[3];
sx q[3];
rz(-0.20488258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4591111) q[0];
sx q[0];
rz(-0.24523188) q[0];
sx q[0];
rz(0.58897513) q[0];
rz(-2.466195) q[1];
sx q[1];
rz(-2.1928936) q[1];
sx q[1];
rz(-1.4640456) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9927514) q[0];
sx q[0];
rz(-2.0939079) q[0];
sx q[0];
rz(-0.4776202) q[0];
x q[1];
rz(1.6590933) q[2];
sx q[2];
rz(-0.73511926) q[2];
sx q[2];
rz(0.59150782) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3551819) q[1];
sx q[1];
rz(-2.1568598) q[1];
sx q[1];
rz(1.5152452) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6028529) q[3];
sx q[3];
rz(-2.3711088) q[3];
sx q[3];
rz(-0.24365261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7318763) q[2];
sx q[2];
rz(-1.7475374) q[2];
sx q[2];
rz(-2.0173006) q[2];
rz(-0.8479979) q[3];
sx q[3];
rz(-0.8046937) q[3];
sx q[3];
rz(-3.1112352) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59415862) q[0];
sx q[0];
rz(-1.1347329) q[0];
sx q[0];
rz(-1.5190079) q[0];
rz(2.2183954) q[1];
sx q[1];
rz(-2.1122439) q[1];
sx q[1];
rz(-1.5069638) q[1];
rz(0.53096622) q[2];
sx q[2];
rz(-0.35561564) q[2];
sx q[2];
rz(1.8993062) q[2];
rz(0.93529978) q[3];
sx q[3];
rz(-2.7705396) q[3];
sx q[3];
rz(-2.0075575) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
