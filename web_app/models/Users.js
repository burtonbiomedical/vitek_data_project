const mongoose = require('mongoose');
const Schema = mongoose.Schema;

const UserSchema = new Schema({
  //We create a property called title, which is mapped to an object that is a string
  //it is also a requirement, so it cannot be left blank
  //and you do this for each property of your schema
  name:{
    type: String,
    required: true
  },
  email:{
    type: String,
    required: true
  },
  admin: {
    type: Boolean,
    required: true
  },
  password:{
    type: String,
    required: true
  },
  date: {
    type: Date,
    default: Date.now
  }
});
//Create schema with name 'users'
mongoose.model('users', UserSchema);