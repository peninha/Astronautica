clearscreen.
print "######### JACY 3 #########" at (0,5).
print "Preparando arquivos..." at (0,6).
wait 2.
print "Copiando 'roll_burn.ks' para sonda..." at (0,8).
copypath("0:/jacy4/roll_burn.ks", "sonda:/").
wait 2.
print "Copiando 'roll.ks' para o segundo estágio..." at (0,9).
copypath("0:/jacy4/roll.ks", "snd_stage:/").
wait 3.
print "Fim da preparação!" at (0,11).
wait 1.
print "Boa sorte com o impacto lunar!" at (0,12).
wait 1.
print "Vai FILHÃÃÃOOO!" at (0,13).
